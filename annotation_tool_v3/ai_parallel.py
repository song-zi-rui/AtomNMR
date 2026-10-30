"""Bounded molecule concurrency and versioned official DeepSeek tariff windows."""
from datetime import datetime, timezone, timedelta
from urllib.parse import urlsplit

DEFAULTS = {'parallelism': 1, 'parallel_strategy': 'fixed'}
MAX_PARALLELISM = 64
MAX_WORKERS = 64
PRICING_SOURCE = 'https://api-docs.deepseek.com/quick_start/pricing/'
PRICING_CHECKED = '2026-10-29'
BEIJING = timezone(timedelta(hours=8))
# State Council holiday notice, 2025-11-04. Weekends remain off-peak even on
# adjusted working days, per the provider's Monday-Friday tariff definition.
HOLIDAYS = {2026: [('01-01', '01-03'), ('02-15', '02-23'), ('04-04', '04-06'),
                   ('05-01', '05-05'), ('06-19', '06-21'), ('09-25', '09-27'),
                   ('10-01', '10-07')]}


def validate(value):
    n = value.get('parallelism', 1)
    strategy = value.get('parallel_strategy', 'fixed')
    if type(n) is not int or not 1 <= n <= MAX_PARALLELISM:
        raise ValueError('并行数量应为 1–64 的整数')
    if strategy not in ('fixed', 'deepseek_offpeak', 'deepseek_offpeak_only'):
        raise ValueError('并行策略无效')
    return {'parallelism': n, 'parallel_strategy': strategy}


def official(config):
    parsed = urlsplit(config.get('base_url', ''))
    return (parsed.scheme == 'https' and parsed.hostname == 'api.deepseek.com'
            and parsed.port in (None, 443) and config.get('api_adapter') == 'deepseek'
            and config.get('model') in ('deepseek-flash', 'deepseek-v4-pro',
                                       'deepseek-v4-flash', 'deepseek-v4-flash-vision-exp'))


def capacity(config, now=None):
    settings = validate(config)
    current = (now or datetime.now(timezone.utc)).astimezone(BEIJING)
    holiday = any(a <= current.strftime('%m-%d') <= b for a, b in HOLIDAYS.get(current.year, []))
    peak = current.weekday() < 5 and not holiday and (9 <= current.hour < 12 or 14 <= current.hour < 18)
    enabled = settings['parallel_strategy'] in ('deepseek_offpeak', 'deepseek_offpeak_only') and official(config)
    n = settings['parallelism']
    only = enabled and settings['parallel_strategy'] == 'deepseek_offpeak_only'
    effective = (0 if peak else n) if only else (16 if enabled and not peak else n)
    return {**settings, 'effective_parallelism': effective,
            'schedule_active': enabled, 'tariff': ('peak' if peak else 'offpeak') if enabled else None,
            'timezone': 'Asia/Shanghai', 'pricing_checked': PRICING_CHECKED,
            'calendar_complete': current.year in HOLIDAYS,
            'note': (('仅谷价运行：峰价不发新请求，谷价按设定数量并行。' if only else '谷价固定 16 并行；峰价时恢复设定数量。')+'切换时已发出的请求继续完成。'
                     + ('' if current.year in HOLIDAYS else '本年度节假日日历尚未更新，工作日按普通工作日保守处理。')) if enabled else
                    '固定并行；谷价策略仅适用于已核实的 DeepSeek 官方端点和模型。'}
