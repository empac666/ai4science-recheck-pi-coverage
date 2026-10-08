import hashlib
import json
from pathlib import Path
import numpy as np


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def clean(x):
    if isinstance(x,dict): return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [clean(v) for v in x]
    if isinstance(x,np.ndarray): return clean(x.tolist())
    if isinstance(x,(np.bool_,)): return bool(x)
    if isinstance(x,np.integer): return int(x)
    if isinstance(x,(float,np.floating)):
        if np.isnan(x): return None
        if np.isposinf(x): return 'Infinity'
        if np.isneginf(x): return '-Infinity'
        return float(x)
    return x


def write_json(path,data):
    Path(path).write_text(json.dumps(clean(data),ensure_ascii=False,indent=2,allow_nan=False)+'\n')
