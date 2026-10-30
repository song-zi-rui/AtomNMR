"""Per-origin CA trust with an explicit legacy-format exception, never verify=False."""
import hashlib
import ssl
from pathlib import Path
from urllib.parse import urlsplit

TLS_DEFAULTS={'tls_policy':'system','ca_file':'','ca_sha256':'','tls_origin':''}


def origin(url):
    parsed=urlsplit(url)
    if parsed.scheme!='https' or not parsed.hostname or parsed.username or parsed.password:
        raise ValueError('指定CA的接口必须使用HTTPS地址')
    return 'https://'+parsed.hostname.lower()+':'+str(parsed.port or 443)


def certificate_fingerprint(path):
    target=Path(path)
    if not target.is_absolute() or not target.is_file():raise ValueError('CA证书文件不存在，请填写绝对路径')
    if target.stat().st_size>128*1024:raise ValueError('CA证书文件过大')
    try:
        text=target.read_text(encoding='utf-8-sig')
        if text.count('-----BEGIN CERTIFICATE-----')!=1 or 'PRIVATE KEY' in text:
            raise ValueError('请选择单张公开PEM格式CA证书，不接受私钥或混合文件')
        der=ssl.PEM_cert_to_DER_cert(text.strip())
        ssl.create_default_context(cafile=str(target))  # Parse the trust file without making a connection.
    except (OSError,UnicodeError,ssl.SSLError,ValueError) as exc:
        raise ValueError('CA文件不是有效的公开PEM证书：'+str(exc)) from None
    return hashlib.sha256(der).hexdigest()


def validate_config(config):
    policy=config.get('tls_policy','system')
    if policy not in ('system','custom_ca','legacy_ca'):raise ValueError('TLS证书策略无效')
    if policy=='system':return dict(TLS_DEFAULTS)
    bound=origin(config['base_url']);saved=config.get('tls_origin','')
    if saved and saved!=bound:raise ValueError('此CA配置已绑定其他接口地址；更换服务器前请先切回系统默认验证，再明确配置新CA')
    path=config.get('ca_file','')
    if not isinstance(path,str):raise ValueError('CA文件路径无效')
    fingerprint=certificate_fingerprint(path)
    expected=config.get('ca_sha256','')
    if expected and expected!=fingerprint:raise ValueError('CA证书内容与已确认指纹不一致；请重新核实证书后配置')
    return {'tls_policy':policy,'ca_file':str(Path(path).resolve()),'ca_sha256':fingerprint,'tls_origin':bound}


def context_for(base_url,settings):
    if not settings or settings.get('tls_policy','system')=='system':return None
    if not settings.get('tls_origin') or not settings.get('ca_sha256'):
        raise ValueError('自定义CA尚未绑定接口与证书指纹，请先保存API配置')
    cfg=validate_config({**settings,'base_url':base_url})
    context=ssl.create_default_context(cafile=cfg['ca_file'])
    if cfg['tls_policy']=='legacy_ca':context.verify_flags &= ~ssl.VERIFY_X509_STRICT
    assert context.verify_mode==ssl.CERT_REQUIRED and context.check_hostname
    return context


def request_options(config):
    if config.get('tls_policy','system')=='system':return {}
    return {'tls_config':{k:config.get(k,v) for k,v in TLS_DEFAULTS.items()}}
