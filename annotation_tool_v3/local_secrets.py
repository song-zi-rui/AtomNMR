"""Protect local API credentials with Windows DPAPI on this computer."""
from __future__ import annotations

import base64
import ctypes
import os
from ctypes import wintypes


MACHINE_PREFIX = 'dpapi-machine-v1:'
UI_FORBIDDEN = 0x1
LOCAL_MACHINE = 0x4


class Blob(ctypes.Structure):
    _fields_ = [('size', wintypes.DWORD), ('data', ctypes.POINTER(ctypes.c_ubyte))]


def _crypt(data, decrypt=False, *, machine=False):
    if os.name != 'nt':
        raise ValueError('本机密钥保存目前需要 Windows DPAPI')
    source_buffer = ctypes.create_string_buffer(data)
    source = Blob(len(data), ctypes.cast(source_buffer, ctypes.POINTER(ctypes.c_ubyte)))
    destination = Blob()
    crypt32 = ctypes.WinDLL('crypt32', use_last_error=True)
    kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)
    kernel32.LocalFree.argtypes = [ctypes.c_void_p]
    kernel32.LocalFree.restype = ctypes.c_void_p
    fn = crypt32.CryptUnprotectData if decrypt else crypt32.CryptProtectData
    second_type = ctypes.POINTER(wintypes.LPWSTR) if decrypt else wintypes.LPCWSTR
    fn.argtypes = [ctypes.POINTER(Blob), second_type, ctypes.POINTER(Blob), ctypes.c_void_p,
                   ctypes.c_void_p, wintypes.DWORD, ctypes.POINTER(Blob)]
    fn.restype = wintypes.BOOL
    try:
        flags = UI_FORBIDDEN | (LOCAL_MACHINE if machine and not decrypt else 0)
        if not fn(ctypes.byref(source), None if decrypt else 'AtomNMR DeepSeek API Key',
                  None, None, None, flags, ctypes.byref(destination)):
            action = '读取' if decrypt else '保存'
            raise ValueError('Windows 无法' + action + '本机加密密钥')
        return ctypes.string_at(destination.data, destination.size)
    finally:
        if not decrypt:
            ctypes.memset(source.data, 0, source.size)
        if destination.data:
            if decrypt:
                ctypes.memset(destination.data, 0, destination.size)
            kernel32.LocalFree(ctypes.cast(destination.data,ctypes.c_void_p))


def protect(value):
    encrypted = _crypt(value.encode('utf-8'), machine=True)
    return MACHINE_PREFIX + base64.b64encode(encrypted).decode('ascii')


def unprotect(value):
    try:
        encrypted = base64.b64decode(value.removeprefix(MACHINE_PREFIX),validate=True)
        return _crypt(encrypted,decrypt=True).decode('utf-8')
    except (UnicodeError, ValueError) as exc:
        if value.startswith(MACHINE_PREFIX):
            raise ValueError('本机密钥无法解密。请在保存它的电脑上使用，或重新输入一次密钥') from exc
        raise ValueError('旧密钥无法解密。请用保存它的 Windows 账号启动一次以迁移，或重新输入一次密钥') from exc
