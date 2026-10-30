"""Estimate production latency from comparable, completed official calls only."""
import json
import math
from . import ai_parallel


def estimate(c,config):
    successes=[];timeouts=[]
    for row in c.execute('''SELECT i.elapsed,i.state,i.request_config_json,i.diagnostics_json,j.config_json
        FROM auto_items i JOIN auto_jobs j ON j.id=i.job_id
        WHERE i.elapsed IS NOT NULL ORDER BY i.updated_at DESC,i.id DESC'''):
        request=json.loads(row['request_config_json'] or row['config_json'])
        if request.get('local_response_reuse') or not ai_parallel.official(request):continue
        if any(request.get(k)!=config.get(k) for k in ('model','thinking','reasoning_effort','wire_format')):continue
        if row['state']=='done':successes.append(float(row['elapsed']))
        elif json.loads(row['diagnostics_json'] or '{}').get('error_kind')=='total_timeout':timeouts.append(float(row['elapsed']))
    def percentile(values,q):
        if not values:return None
        values=sorted(values);index=(len(values)-1)*q;lo=math.floor(index);hi=math.ceil(index)
        return values[lo]+(values[hi]-values[lo])*(index-lo)
    # Recent tail captures current serving speed; the client timeout is an
    # additional conservative floor, not a claim that non-streaming is bounded.
    recent=successes[:500]
    upper=max([300,float(config.get('timeout',600)),*recent[:500],*timeouts[:50]])
    statistical=1.25*(percentile(recent,.99) or 0)
    guard=math.ceil(max(upper,statistical)+60)
    return {'samples':len(successes),'mean_seconds':sum(successes)/len(successes) if successes else None,
        'median_seconds':percentile(successes,.5),'p95_seconds':percentile(successes,.95),
        'p99_seconds':percentile(successes,.99),'max_seconds':max(successes) if successes else None,
        'recent_samples':len(recent),'guard_seconds':guard,'margin_seconds':60,
        'method':'max(最近500次最慢耗时, 最近P99×1.25, 配置等待上限, 近期总超时耗时, 300秒)+60秒；不含本地复用响应',
        'guaranteed':False}
