"""Recognize one missing OUTERMOST object brace; never invent answer fields."""
import json
from .ai_protocol import unique_object_pairs


def close_missing_root_object(text):
    if not isinstance(text,str):return None
    stripped=text.strip()
    if not stripped.startswith('{') or not stripped.endswith(('}',']')):return None
    stack=[];quoted=False;escaped=False
    for char in stripped:
        if quoted:
            if escaped:escaped=False
            elif char=='\\':escaped=True
            elif char=='"':quoted=False
            continue
        if char=='"':quoted=True
        elif char in '{[':stack.append('}' if char=='{' else ']')
        elif char in '}]':
            if not stack or stack.pop()!=char:return None
    if quoted or stack!=['}']:return None
    candidate=stripped+'}'
    try:json.loads(candidate,object_pairs_hook=unique_object_pairs)
    except (ValueError,TypeError):return None
    return candidate
