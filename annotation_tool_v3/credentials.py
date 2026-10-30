"""Named, locally encrypted API credentials with explicit lifecycle actions."""
import hashlib
import json
import uuid
from contextlib import closing

from . import core as db, local_secrets


def initialize(c):
    c.executescript('''
    CREATE TABLE IF NOT EXISTS api_credentials(
      id TEXT PRIMARY KEY,name TEXT NOT NULL,provider TEXT NOT NULL,encrypted TEXT NOT NULL,
      fingerprint TEXT NOT NULL,suffix TEXT NOT NULL,enabled INTEGER NOT NULL,
      created_at TEXT NOT NULL,updated_at TEXT NOT NULL,last_used_at TEXT,
      verification TEXT,verified_at TEXT);
    CREATE UNIQUE INDEX IF NOT EXISTS credential_fingerprint ON api_credentials(provider,fingerprint);
    ''')
    old=c.execute("SELECT value FROM settings WHERE key='deepseek_api_key_dpapi'").fetchone()
    if old:
        # Preserve a legacy ciphertext even when it belongs to another Windows user.
        try:
            secret=local_secrets.unprotect(old[0]);fingerprint=hashlib.sha256(secret.encode()).hexdigest();suffix=secret[-4:]
        except ValueError:
            fingerprint='legacy:'+hashlib.sha256(old[0].encode()).hexdigest();suffix='待更新'
        identifier=uuid.uuid4().hex
        c.execute('INSERT OR IGNORE INTO api_credentials VALUES(?,?,?,?,?,?,1,?,?,NULL,NULL,NULL)',
                  (identifier,'原有 DeepSeek 密钥','deepseek',old[0],fingerprint,suffix,db.now(),db.now()))
        actual=c.execute("SELECT id FROM api_credentials WHERE provider='deepseek' AND fingerprint=?",(fingerprint,)).fetchone()[0]
        c.execute("INSERT OR REPLACE INTO settings VALUES('active_credential',?)",(json.dumps(actual),))
        c.execute("DELETE FROM settings WHERE key='deepseek_api_key_dpapi'")
    c.commit()


class Vault:
    def __init__(self,path,connector):
        self.path,self.connect=path,connector

    def active(self,c):
        r=c.execute("SELECT value FROM settings WHERE key='active_credential'").fetchone()
        return json.loads(r[0]) if r else None

    def status(self):
        with closing(self.connect(self.path)) as c:
            active=self.active(c)
            rows=c.execute('SELECT id,name,provider,suffix,enabled,created_at,updated_at,last_used_at,verification,verified_at FROM api_credentials ORDER BY created_at,id').fetchall()
        return {'active_id':active,'credentials':[{**dict(r),'masked_key':'•••• '+r['suffix'],'active':r['id']==active} for r in rows]}

    def read(self,identifier=None):
        with closing(self.connect(self.path)) as c:
            r=c.execute('SELECT * FROM api_credentials WHERE id=?',(identifier or self.active(c),)).fetchone()
        if not r:return '',None,None
        if not r['enabled']:return '',dict(r),'当前密钥已停用，请切换或启用'
        try:return local_secrets.unprotect(r['encrypted']),dict(r),None
        except ValueError as exc:return '',dict(r),str(exc)

    def change(self,payload):
        action=payload.get('action')
        with closing(self.connect(self.path)) as c:
            c.execute('BEGIN IMMEDIATE')
            active=self.active(c)
            identifier=payload.get('id') or active
            row=c.execute('SELECT * FROM api_credentials WHERE id=?',(identifier,)).fetchone()
            if action in ('create','replace'):
                provider=payload.get('provider') or (row['provider'] if row else 'deepseek')
                if provider not in ('deepseek','qwen','openai-compatible'):raise ValueError('密钥供应方无效')
                secret=payload.get('api_key')
                if not isinstance(secret,str) or not 8<=len(secret)<=512 or any(x.isspace() or ord(x)<32 for x in secret):
                    raise ValueError('请填写有效密钥，不要使用脱敏显示值')
                name=str(payload.get('name') or (row['name'] if row else '')).strip()
                if not name or len(name)>80:raise ValueError('请填写 1–80 字的密钥名称')
                if action=='replace' and not row:raise ValueError('要替换的密钥不存在')
                fingerprint=hashlib.sha256(secret.encode()).hexdigest()
                duplicate=c.execute('SELECT id FROM api_credentials WHERE provider=? AND fingerprint=?',(provider,fingerprint)).fetchone()
                if duplicate and (action=='create' or duplicate[0]!=identifier):raise ValueError('此密钥已保存，请直接选择它')
                encrypted=local_secrets.protect(secret)
                if action=='create':
                    identifier=uuid.uuid4().hex
                    c.execute('INSERT INTO api_credentials VALUES(?,?,?,?,?,?,1,?,?,NULL,NULL,NULL)',
                              (identifier,name,provider,encrypted,fingerprint,secret[-4:],db.now(),db.now()))
                else:
                    c.execute('UPDATE api_credentials SET name=?,provider=?,encrypted=?,fingerprint=?,suffix=?,enabled=1,updated_at=?,verification=NULL,verified_at=NULL WHERE id=?',
                              (name,provider,encrypted,fingerprint,secret[-4:],db.now(),identifier))
                if payload.get('activate',True) or active==identifier:
                    c.execute("INSERT OR REPLACE INTO settings VALUES('active_credential',?)",(json.dumps(identifier),))
            elif action=='activate':
                if not row or not row['enabled']:raise ValueError('请先启用该密钥')
                local_secrets.unprotect(row['encrypted'])
                c.execute("INSERT OR REPLACE INTO settings VALUES('active_credential',?)",(json.dumps(identifier),))
            elif action in ('disable','enable','remove','rename'):
                if not row:raise ValueError('密钥不存在')
                if action=='remove':
                    c.execute('DELETE FROM api_credentials WHERE id=?',(identifier,))
                    if active==identifier:c.execute("DELETE FROM settings WHERE key='active_credential'")
                elif action=='rename':
                    name=str(payload.get('name') or '').strip()
                    if not name or len(name)>80:raise ValueError('密钥名称应为 1–80 字')
                    c.execute('UPDATE api_credentials SET name=?,updated_at=? WHERE id=?',(name,db.now(),identifier))
                else:c.execute('UPDATE api_credentials SET enabled=?,updated_at=? WHERE id=?',(int(action=='enable'),db.now(),identifier))
            else:raise ValueError('未知密钥管理操作')
            c.commit()
        return self.status()

    def mark(self,identifier,*,verification=None,used=False):
        if not identifier:return
        with closing(self.connect(self.path)) as c:
            if used:c.execute('UPDATE api_credentials SET last_used_at=? WHERE id=?',(db.now(),identifier))
            if verification is not None:c.execute('UPDATE api_credentials SET verification=?,verified_at=? WHERE id=?',(verification,db.now(),identifier))
            c.commit()
