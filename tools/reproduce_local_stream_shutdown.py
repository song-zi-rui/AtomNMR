"""Loopback-only TLS control: distinguish local deadline from server failure."""
import json
from pathlib import Path
import socket
import ssl
import subprocess
import tempfile
import threading
import time
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from urllib.request import build_opener,HTTPSHandler,ProxyHandler

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/transport_root_cause_20261031'
OUT.mkdir(exist_ok=True)

def main():
    openssl=Path('C:/Users/ray/software/miniconda3/Library/bin/openssl.exe')
    results=[]
    with tempfile.TemporaryDirectory(prefix='local-tls-',dir=OUT) as directory:
        folder=Path(directory);cert=folder/'loopback.crt';key=folder/'loopback.key'
        config=folder/'openssl.cnf'
        config.write_text('[req]\ndistinguished_name=dn\nx509_extensions=v3\n[dn]\nCN=localhost\n[v3]\nbasicConstraints=critical,CA:TRUE\nkeyUsage=critical,digitalSignature,keyEncipherment,keyCertSign\nsubjectKeyIdentifier=hash\nauthorityKeyIdentifier=keyid:always\n',encoding='ascii')
        generated=subprocess.run([str(openssl),'req','-config',str(config),'-x509','-newkey','rsa:2048','-nodes','-days','1',
            '-subj','/CN=localhost','-addext','subjectAltName=IP:127.0.0.1',
            '-keyout',str(key),'-out',str(cert)],stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        if generated.returncode:raise ValueError('Local certificate generation failed: '+generated.stderr.decode(errors='replace')[:1500])
        for name in ['normal_completion','local_deadline_during_stall','local_shutdown_before_next_read','remote_eof_without_finish']:
            finished=threading.Event()
            class Handler(BaseHTTPRequestHandler):
                protocol_version='HTTP/1.1'
                def log_message(self,*args):pass
                def do_GET(self):
                    self.send_response(200);self.send_header('Content-Type','text/event-stream');self.send_header('Connection','close');self.end_headers()
                    self.wfile.write(b'data: {"choices":[{"delta":{"content":"OK"},"finish_reason":null}]}\n\n');self.wfile.flush()
                    if name=='normal_completion':
                        time.sleep(.05);self.wfile.write(b'data: {"choices":[{"delta":{},"finish_reason":"stop"}]}\n\ndata: [DONE]\n\n');self.wfile.flush()
                    elif name.startswith('local_'):finished.wait(7)
                    self.close_connection=True
            server=ThreadingHTTPServer(('127.0.0.1',0),Handler);server.daemon_threads=True
            context=ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER);context.load_cert_chain(cert,key)
            server.socket=context.wrap_socket(server.socket,server_side=True)
            thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
            client_context=ssl.create_default_context(cafile=str(cert))
            opener=build_opener(ProxyHandler({}),HTTPSHandler(context=client_context))
            started=time.perf_counter();response=None;trace=[];flag=threading.Event()
            def tick(label):trace.append({'elapsed_s':round(time.perf_counter()-started,6),'event':label})
            def deadline():
                tick('local_deadline_reached');flag.set()
                if response is not None:
                    try:
                        tick('local_shutdown_called')
                        response.fp.raw._sock.shutdown(socket.SHUT_RDWR)
                    except (OSError,AttributeError) as exc:tick('shutdown_'+type(exc).__name__)
            timer=threading.Timer(1.5,deadline);timer.daemon=True
            result={'case':name,'tls_certificate_verified':True,'external_requests':0}
            timer.start()
            try:
                response=opener.open('https://127.0.0.1:'+str(server.server_port),timeout=5);tick('http_200')
                while True:
                    line=response.readline()
                    if not line:tick('eof');result['outcome']='eof';break
                    if b'"content":"OK"' in line:
                        tick('answer_received')
                        if name=='local_shutdown_before_next_read':
                            flag.wait(2);tick('read_again_after_shutdown')
                    if b'"finish_reason":"stop"' in line:tick('finish_received')
                    if b'[DONE]' in line:tick('done_received');result['outcome']='complete';break
            except Exception as exc:
                tick('exception_received');result.update(outcome='exception',exception_type=type(exc).__name__,winerror=getattr(exc,'winerror',None),message=str(exc))
            finally:
                timer.cancel();finished.set()
                if response is not None:response.close()
                server.shutdown();server.server_close();thread.join(2)
            result.update(deadline_fired=flag.is_set(),trace=trace);results.append(result)
    (OUT/'local_tls_controls.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(results,ensure_ascii=True))

if __name__=='__main__':main()
