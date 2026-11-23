"""Classify failures for unattended production annotation retries."""
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from typing import NamedTuple

from . import ai_annotation as ai


CONNECTION_DELAYS = (30, 90)
UNCERTAIN_DELAYS = (60,)
OUTPUT_DELAYS = (5,)
STREAM_ERROR_DELAYS = (30,)
LEGACY_STREAM_ERROR = '流式接口报告服务错误；本批次停止，不自动重试。'
TRANSIENT_HTTP = {408, 425, 429, 500, 502, 503, 504, 520, 521, 522, 523, 524}
TRANSIENT_NETWORK = {'dns', 'connection_refused', 'tls_protocol', 'network', 'timeout',
                     'total_timeout', 'connection_reset', 'remote_disconnected', 'incomplete_read', 'incomplete_stream'}
QUOTA_ERRORS = ('insufficient_quota', 'insufficient_balance', 'billing_hard_limit',
                'budget_exceeded', 'budgetexceedederror', '余额不足', '欠费', 'quota exhausted')


class Retry(NamedTuple):
    reason: str
    delays: tuple
    cool_down: bool
    pause_on_exhaustion: bool


def classify(exc):
    diagnostics = getattr(exc, 'diagnostics', {}) or {}
    kind = diagnostics.get('error_kind')
    status = diagnostics.get('provider_status') if kind == 'stream_service_error' else diagnostics.get('http_status')
    uncertain = diagnostics.get('request_outcome_uncertain') is True
    info = diagnostics.get('provider_error') or {}
    desc = ' '.join(str(v) for v in info.values()).lower()
    if any(reason in desc for reason in QUOTA_ERRORS):
        return None
    if status in TRANSIENT_HTTP:
        delays = UNCERTAIN_DELAYS if uncertain else CONNECTION_DELAYS
        label = 'SSE' if kind == 'stream_service_error' else 'HTTP'
        return Retry(f'{label} {status}', delays, True, status in (429, 503, 521, 523))
    if kind == 'stream_service_error' and status is None:
        stream = diagnostics.get('stream_transport') or {}
        delays = UNCERTAIN_DELAYS if stream.get('reasoning_characters') or stream.get('content_characters') else STREAM_ERROR_DELAYS
        return Retry(kind, delays, True, False)
    if kind in TRANSIENT_NETWORK:
        delays = UNCERTAIN_DELAYS if uncertain else CONNECTION_DELAYS
        return Retry(kind, delays, True, not uncertain)
    if kind == 'repetitive_output':
        return Retry(kind, UNCERTAIN_DELAYS if uncertain else OUTPUT_DELAYS, False, False)
    if isinstance(exc, ai.AIResponseError):
        if diagnostics.get('finish_reason') in (None, 'stop'):
            return Retry(kind or 'invalid_output', OUTPUT_DELAYS, False, False)
    return None


def retry_after_seconds(header):
    if not isinstance(header, str):
        return 0
    if header.isascii() and header.isdigit():
        return int(header)
    try:
        date = parsedate_to_datetime(header)
        if date.tzinfo is None:
            date = date.replace(tzinfo=timezone.utc)
        return max(0, (date - datetime.now(timezone.utc)).total_seconds())
    except (ValueError, TypeError, OverflowError):
        return 0
