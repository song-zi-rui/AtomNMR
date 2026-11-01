"""Interrupt local HTTP waits; does not claim to cancel remote model execution."""
import http.client
import ipaddress
import socket
import threading
import time
from contextlib import contextmanager
from contextvars import ContextVar
from functools import partial
from urllib.request import HTTPHandler,HTTPSHandler

_current=ContextVar('ai_request_control',default=None)


class APIRequestCancelled(ValueError):
    def __init__(self):
        super().__init__('用户暂停：已停止本地等待；远端可能仍在处理，结果未确认，不自动重试此分子。')
        self.usage={}
        self.uncertain=True
        self.diagnostics={'user_cancelled':True,'request_outcome_uncertain':True,'transport_error':False}


class RequestControl:
    def __init__(self):
        self._lock=threading.Lock();self._sockets=set();self._cancelled=False;self._finished=False
        self._phase=None;self._application_send_started=False

    def stage(self,phase,application_send_started=False):
        with self._lock:
            self._phase=phase
            self._application_send_started|=application_send_started

    def transport_state(self):
        with self._lock:
            return {'failure_stage':self._phase,'application_send_started':self._application_send_started}

    def check(self):
        with self._lock:cancelled=self._cancelled and not self._finished
        if cancelled:raise APIRequestCancelled()

    def register(self,sock):
        with self._lock:
            cancelled=self._cancelled and not self._finished
            if not cancelled:self._sockets.add(sock)
        if cancelled:
            self._close(sock);raise APIRequestCancelled()

    def unregister(self,sock):
        with self._lock:self._sockets.discard(sock)

    @staticmethod
    def _close(sock):
        try:sock.shutdown(socket.SHUT_RDWR)
        except (OSError,ValueError):pass
        # socket.makefile() retains IO references. On Windows ordinary close()
        # defers closesocket(), so a blocked read can otherwise survive pause.
        try:
            handle=sock.detach()
            if handle!=-1:socket.close(handle)
        except (OSError,ValueError):pass

    def cancel(self):
        with self._lock:self._cancelled=True;sockets=tuple(self._sockets)
        for sock in sockets:self._close(sock)

    def seal(self):
        # A fully received response wins a later pause and may still be saved.
        with self._lock:
            if self._cancelled and not self._finished:raise APIRequestCancelled()
            self._finished=True

    def close(self):
        with self._lock:sockets=tuple(self._sockets);self._sockets.clear()
        for sock in sockets:self._close(sock)


@contextmanager
def capture_request_control(control):
    token=_current.set(control)
    try:yield control
    finally:_current.reset(token)


def current_request_control():return _current.get()


def _addresses(host,port,control,deadline):
    control.stage('dns')
    try:ipaddress.ip_address(host)
    except ValueError:
        done=threading.Event();result=[];errors=[]
        def lookup():
            try:result.extend(socket.getaddrinfo(host,port,0,socket.SOCK_STREAM))
            except Exception as exc:errors.append(exc)
            finally:done.set()
        threading.Thread(target=lookup,daemon=True,name='ai-dns').start()
        while not done.wait(.05):
            control.check()
            if deadline is not None and time.monotonic()>=deadline:raise TimeoutError('DNS lookup timed out')
        control.check()
        if errors:raise errors[0]
        return result
    return socket.getaddrinfo(host,port,0,socket.SOCK_STREAM,0,socket.AI_NUMERICHOST)


def _connect(address,timeout,source_address,control):
    actual=socket.getdefaulttimeout() if timeout is socket._GLOBAL_DEFAULT_TIMEOUT else timeout
    deadline=time.monotonic()+actual if actual is not None else None
    errors=[]
    for family,kind,protocol,_,target in _addresses(*address,control,deadline):
        control.check();sock=socket.socket(family,kind,protocol);control.register(sock)
        try:
            if deadline is not None:
                left=deadline-time.monotonic()
                if left<=0:raise TimeoutError('connection timed out')
                sock.settimeout(left)
            if source_address:sock.bind(source_address)
            control.stage('tcp_connect')
            sock.connect(target);control.check();return sock
        except Exception as exc:
            control.unregister(sock);RequestControl._close(sock);control.check();errors.append(exc)
    if errors:raise errors[-1]
    raise OSError('no connection address available')


class _ControlledConnection:
    def __init__(self,*args,control,**kwargs):
        self._control=control;self._connecting=False;super().__init__(*args,**kwargs)
        self._create_connection=lambda address,timeout,source_address: _connect(address,timeout,source_address,control)

    def connect(self):
        self._connecting=True
        try:return super().connect()
        finally:self._connecting=False

    def _tunnel(self):
        self._control.stage('proxy_tunnel')
        return super()._tunnel()

    def send(self,data):
        self._control.check()
        if self.sock is None:
            if not self.auto_open:raise http.client.NotConnected()
            self.connect()
        # CONNECT is proxy negotiation, not submission of the annotation payload.
        if not self._connecting:self._control.stage('request_write',application_send_started=True)
        return super().send(data)

    def getresponse(self):
        self._control.stage('response_headers')
        return super().getresponse()


class _HTTPConnection(_ControlledConnection,http.client.HTTPConnection):pass


class _HTTPSConnection(_ControlledConnection,http.client.HTTPSConnection):
    def connect(self):
        self._connecting=True
        try:
            http.client.HTTPConnection.connect(self)
            self._control.check();self._control.stage('tls_handshake')
            raw=self.sock
            self.sock=self._context.wrap_socket(raw,server_hostname=self._tunnel_host or self.host,do_handshake_on_connect=False)
            self._control.unregister(raw);self._control.register(self.sock)
            self.sock.do_handshake();self._control.check()
        finally:self._connecting=False


class ControlledHTTPHandler(HTTPHandler):
    def __init__(self,control):super().__init__();self.control=control
    def http_open(self,req):return self.do_open(partial(_HTTPConnection,control=self.control),req)


class ControlledHTTPSHandler(HTTPSHandler):
    def __init__(self,control,context=None):super().__init__(context=context);self.control=control
    def https_open(self,req):return self.do_open(partial(_HTTPSConnection,control=self.control),req,context=self._context)
