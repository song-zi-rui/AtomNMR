"""Streaming transport status. Never retains model reasoning text."""
import json
import re
import socket
import ssl
import time
from contextlib import contextmanager
from contextvars import ContextVar
from email.utils import parsedate_to_datetime
from .ai_cancellation import current_request_control

_observer = ContextVar('ai_transport_progress', default=None)


class APIRequestError(ValueError):
    def __init__(self, message, *, status=None, uncertain=False, usage=None):
        super().__init__(message)
        self.status = status
        self.uncertain = uncertain
        self.should_pause_batch = True
        self.usage = usage or {}
        self.diagnostics = {'http_status': status, 'request_outcome_uncertain': uncertain,
                            'transport_error': True}


def _error_text(value, key='', limit=500):
    """Only small, selected error fields are retained; never the full response."""
    if not isinstance(value,(str,int,float)) or isinstance(value,bool):return ''
    text=str(value)
    if key:text=text.replace(key,'[REDACTED]')
    text=re.sub(r'(?i)\bBearer\s+[^\s\"\'<>]+','Bearer [REDACTED]',text)
    text=re.sub(r'(?i)\b(?:sk|api[_-]?key)[_-][A-Za-z0-9_-]{12,}','[REDACTED]',text)
    text=re.sub(r'[\x00-\x1f\x7f]+',' ',text)
    return text.strip()[:limit]


def http_error(response, *, key='', elapsed_s=None):
    """Expose provider reasons and retry guidance without resubmitting a request."""
    status=response.code
    headers={str(k).lower():str(v) for k,v in (response.headers or {}).items()}
    body=b''
    try:body=response.read(16385)
    except (OSError,ValueError):pass
    finally:response.close()
    provider={}
    if len(body)<=16384:
        try:
            value=json.loads(body)
            if isinstance(value,dict):
                detail=value.get('error',value.get('detail',value))
                if isinstance(detail,dict):
                    provider={k:_error_text(detail.get(k),key,600 if k=='message' else 100)
                              for k in ('message','code','type','param')}
                elif isinstance(detail,str):provider={'message':_error_text(detail,key,600)}
                if not provider.get('message'):provider['message']=_error_text(value.get('message'),key,600)
                provider['request_id']=_error_text(value.get('request_id'),key,160)
        except (ValueError,UnicodeError):
            # Do not archive HTML, stack traces, echoed prompts or arbitrary bodies.
            if 'text/plain' in headers.get('content-type',''):
                text=body.decode('utf-8',errors='replace').strip()
                if not text.startswith('<'):provider['message']=_error_text(text,key,600)
    messages={401:'API Key 无效',402:'接口余额不足',429:'接口拒绝请求（HTTP 429）',400:'模型或请求参数不受支持',
              404:'API路径不存在（HTTP 404），请核对Base URL是否包含 /v1 等正确路径',
              301:'API地址要求重定向，未转发密钥；请填写最终API地址',302:'API地址要求重定向，未转发密钥；请填写最终API地址',
              524:'网关等待模型响应超时（HTTP 524）；本地增大超时不能延长网关时限。可尝试流式接收或关闭思考后另开测试。',
              502:'模型上游服务异常（HTTP 502），请稍后再测试',503:'模型服务暂不可用或过载（HTTP 503）',
              504:'网关等待上游超时（HTTP 504）'}
    message=messages.get(status,'模型接口错误 HTTP '+str(status))
    if provider.get('message'):message+='；服务端原因：'+provider['message']
    elif status==429:message+='；服务端未提供具体原因，尚不能区分频率、并发或额度限制'
    retry=_error_text(headers.get('retry-after'),key,100)
    if retry:
        if retry.isascii() and retry.isdigit():message+='；服务端建议等待 '+retry+' 秒后再试'
        else:
            try:
                if parsedate_to_datetime(retry):message+='；服务端建议重试时间：'+retry
            except (ValueError,TypeError,OverflowError):retry=''
    request_id=next((_error_text(headers.get(k),key,160) for k in ('x-request-id','request-id','x-correlation-id')
                     if headers.get(k)),provider.get('request_id',''))
    if request_id:message+='；请求编号：'+request_id
    error=APIRequestError(message,status=status,uncertain=status>=500)
    error.diagnostics.update(error_kind='http',provider_error={k:v for k,v in provider.items() if v},
                             request_id=request_id or None,retry_after=retry or None)
    if elapsed_s is not None:error.diagnostics['elapsed_s']=elapsed_s
    rate_headers={k:_error_text(v,key,120) for k,v in headers.items()
                  if k in ('x-ratelimit-limit-requests','x-ratelimit-remaining-requests','x-ratelimit-reset-requests',
                           'x-ratelimit-limit-tokens','x-ratelimit-remaining-tokens','x-ratelimit-reset-tokens')}
    if rate_headers:error.diagnostics['rate_limit_headers']=rate_headers
    return error


def network_error(exc):
    """Keep the failure layer explicit without relaxing TLS or exposing credentials."""
    reason=getattr(exc,'reason',exc)
    if isinstance(reason,ssl.SSLCertVerificationError):
        detail=str(getattr(reason,'verify_message','') or reason).replace('\n',' ')[:240]
        error=APIRequestError('HTTPS证书验证失败，未完成TLS握手；归属请求未提交到模型接口。请检查服务器证书链、CA及访问域名。详情：'+detail,
                              uncertain=False)
        error.diagnostics.update(error_kind='tls_certificate',failure_stage='tls_handshake',
                                 ssl_verify_code=getattr(reason,'verify_code',None),ssl_verify_message=detail)
    elif isinstance(reason,socket.gaierror):
        error=APIRequestError('API主机名无法解析；请求尚未提交。请检查地址、DNS及内网/VPN连接。',uncertain=False)
        error.diagnostics.update(error_kind='dns',failure_stage='dns')
    elif isinstance(reason,ConnectionRefusedError):
        error=APIRequestError('API连接被拒绝；请求尚未提交。请检查服务器端口、服务是否启动及防火墙。',uncertain=False)
        error.diagnostics.update(error_kind='connection_refused',failure_stage='tcp_connect')
    elif isinstance(reason,ssl.SSLError):
        error=APIRequestError('TLS连接异常；请核对HTTP/HTTPS协议及服务器TLS配置。远端处理结果未知，不自动重试。',uncertain=True)
        error.diagnostics.update(error_kind='tls_protocol')
    elif isinstance(reason,TimeoutError):
        error=APIRequestError('API连接或读取超时；远端处理结果未知，不自动重试。请检查服务排队、响应时间和网络。',uncertain=True)
        error.diagnostics.update(error_kind='timeout')
    elif isinstance(reason,ConnectionResetError):
        error=APIRequestError('API连接被中途重置；远端处理结果未知，不自动重试。请检查服务端和代理日志。',uncertain=True)
        error.diagnostics.update(error_kind='connection_reset')
    else:
        error=APIRequestError('连接失败或响应中断；远端处理结果未知，不自动重试，请检查网络与服务日志。',uncertain=True)
        error.diagnostics.update(error_kind='network')
    error.diagnostics['exception_type']=type(reason).__name__
    return error


@contextmanager
def capture_progress(callback):
    token = _observer.set(callback)
    try:
        yield
    finally:
        _observer.reset(token)


def notify_progress(**value):
    callback = _observer.get()
    if callback:
        try:
            callback(value)
        except Exception:
            pass  # A progress-display failure must not interrupt a submitted request.


def should_pause_batch(exc):
    if getattr(exc, 'should_pause_batch', False):
        return True
    # Retain compatibility with provider adapters and old mocks raising ValueError.
    message = str(exc)
    return any(s in message for s in ('HTTP 500', 'HTTP 502', 'HTTP 503', 'HTTP 504', 'HTTP 524',
                                     '连接失败或超时', 'API Key', '余额不足', '限流', '密钥'))


def repetitive_output_period(text):
    """Only flag a short suffix repeated 64 times, not repeated JSON field names."""
    if not isinstance(text,str):return None
    tail=''.join(text[-16384:].split())
    for width in range(3,min(80,len(tail)//64)+1):
        unit=tail[-width:]
        if any(ch.isalnum() for ch in unit) and tail.endswith(unit*64):return width
    return None


def read_sse(response, *, started, timeout):
    content = []
    usage = {}
    meta = {'stream': True, 'chunks': 0, 'content_characters': 0, 'reasoning_characters': 0,
            'ttft_ms': None}
    result = {'choices': [], 'usage': usage}
    finish = None
    done = False
    size = 0
    pending = []
    last_report = 0
    answer_tail = ''
    repeat_checked_at = 0

    def report(force=False):
        nonlocal last_report
        now = time.perf_counter()
        if force or now - last_report >= .75:
            notify_progress(phase='streaming' if meta['ttft_ms'] is not None else 'waiting_first_token', **meta)
            last_report = now

    def event():
        nonlocal finish, done, answer_tail, repeat_checked_at
        if not pending:
            return
        text = '\n'.join(pending)
        pending.clear()
        if text.strip() == '[DONE]':
            done = True
            return
        try:
            value = json.loads(text)
        except (ValueError, UnicodeError):
            raise APIRequestError('流式接口返回了无效数据；本批次停止，请检查网关兼容性。', uncertain=True, usage=usage) from None
        if not isinstance(value, dict):
            raise APIRequestError('流式接口返回格式异常；请检查网关兼容性。', uncertain=True, usage=usage)
        if value.get('error'):
            raise APIRequestError('流式接口报告服务错误；本批次停止，不自动重试。', uncertain=True, usage=usage)
        meta['chunks'] += 1
        for key in ('id', 'model', 'created'):
            if value.get(key) is not None:
                result[key] = value[key]
        if isinstance(value.get('usage'), dict):
            usage.update(value['usage'])
        for choice in value.get('choices') or []:
            if choice.get('index', 0) != 0:
                continue
            delta = choice.get('delta') or {}
            answer, reasoning = delta.get('content'), delta.get('reasoning_content') or delta.get('reasoning')
            if isinstance(answer, str):
                content.append(answer)
                meta['content_characters'] += len(answer)
                answer_tail=(answer_tail+answer)[-16384:]
            if isinstance(reasoning, str):
                meta['reasoning_characters'] += len(reasoning)
            if (answer or reasoning) and meta['ttft_ms'] is None:
                meta['ttft_ms'] = (time.perf_counter() - started) * 1000
            if choice.get('finish_reason') is not None:
                finish = choice['finish_reason']
        if meta['content_characters']>=1024 and meta['content_characters']-repeat_checked_at>=512:
            repeat_checked_at=meta['content_characters']
            period=repetitive_output_period(answer_tail)
            if period:
                report(True)
                error=APIRequestError('模型输出进入连续重复循环；已停止接收并暂停批次，不采用不完整归属。请检查模型或推理参数，不建议仅增加token上限；远端可能仍在运行，不自动重试。',
                                      status=200,uncertain=finish is None,usage=usage)
                error.diagnostics.update(error_kind='repetitive_output',transport_error=False,
                                         repeated_period_characters=period,stream_transport=dict(meta))
                raise error
        report()

    report(True)
    while not done:
        control=current_request_control()
        if control:control.check()
        elapsed=time.perf_counter()-started
        if elapsed > timeout:
            report(True)
            stage=('answer_generation' if meta['content_characters'] else
                   'reasoning' if meta['reasoning_characters'] else 'waiting_first_token')
            error=APIRequestError(
                f"本地总时限 {timeout:g} 秒已到；已接收思考 {meta['reasoning_characters']} 字符、归属答案 {meta['content_characters']} 字符。"
                '已暂停批次，远端是否完成未知，不自动重试。',status=200,uncertain=True,usage=usage)
            error.diagnostics.update(error_kind='total_timeout',failure_stage=stage,transport_error=False,
                                     elapsed_s=elapsed,timeout_s=timeout,stream_transport=dict(meta))
            raise error
        line = response.readline(1024 * 1024 + 1)
        if not line:
            event()
            break
        size += len(line)
        if size > 8 * 1024 * 1024 or len(line) > 1024 * 1024:
            raise APIRequestError('流式响应过大；停止接收并暂停批次。', uncertain=True, usage=usage)
        try:
            text = line.decode('utf-8').rstrip('\r\n')
        except UnicodeError:
            raise APIRequestError('流式响应编码异常；请检查网关兼容性。', uncertain=True, usage=usage) from None
        if not text:
            event()
        elif text.startswith('data:'):
            pending.append(text[5:].lstrip(' '))
        elif text.lstrip().lower().startswith(('<!doctype html', '<html')):
            raise APIRequestError('API地址返回了网页；请核对Base URL的API路径。')
    if finish is None:
        raise APIRequestError('流式连接结束，但未收到完整结束标记；已暂停批次，不自动重复请求。', uncertain=True, usage=usage)
    report(True)
    result.update(choices=[{'index': 0, 'finish_reason': finish, 'message': {'role': 'assistant', 'content': ''.join(content)}}],
                  _transport=meta)
    return result
