"""Local timing statistics; durations are not inferred from token counts."""
import math
import statistics


def timings(values):
    values=sorted(float(v) for v in values if v is not None)
    if not values:return {'n':0,'total_s':0,'mean_s':None,'median_s':None,'p90_s':None,'max_s':None}
    return {'n':len(values),'total_s':sum(values),'mean_s':statistics.mean(values),
            'median_s':statistics.median(values),'p90_s':values[max(0,math.ceil(.9*len(values))-1)],'max_s':values[-1]}
