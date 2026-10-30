"""Protect local API credentials with Windows DPAPI for the current user."""
from __future__ import annotations

import base64
import ctypes
import os
from ctypes import wintypes


class Blob(ctypes.Structure):
    _fields_ = [('size', wintypes.DWORD), ('data', ctypes.POINTER(ctypes.c_ubyte))]


def _crypt(data, decrypt=False):
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
        # UI forbidden; no machine-wide flag, so another Windows user cannot decrypt it.
        if not fn(ctypes.byref(source), None if decrypt else 'AtomNMR DeepSeek API Key',
                  None, None, None, 1, ctypes.byref(destination)):
            action = '读取' if decrypt else '保存'
            raise ValueError('无法' + action + '本机加密密钥；请使用同一 Windows 用户启动标注服务，或重新填写密钥')
        return ctypes.string_at(destination.data, destination.size)
    finally:
        if destination.data:
            kernel32.LocalFree(ctypes.cast(destination.data,ctypes.c_void_p))


def protect(value):
    return base64.b64encode(_crypt(value.encode('utf-8'))).decode('ascii')


def unprotect(value):
    try:
        encrypted = base64.b64decode(value,validate=True)
        return _crypt(encrypted,decrypt=True).decode('utf-8')
    except (UnicodeError, ValueError) as exc:
        raise ValueError('已保存密钥无法解密；请使用原 Windows 用户启动服务或重新输入密钥') from exc
