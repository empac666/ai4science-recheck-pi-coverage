"""Frozen paired bootstrap, compatibility gate, and reporting decision mapping."""
import numpy as np

Q = '在冻结轨道 x* 的观测噪声总体、测试索引 {160,…,799} 且 A 指本稿主约定的范围内，作为对已知原表差距的独立实现检验，'


def bias_diagnostic(point, boot):
    mean = float(np.mean(boot))
    bias = mean-point
    se = float(np.std(boot, ddof=1))
    b = bias/se if se else ('N/A (degenerate)' if bias == 0 else 'ERROR (degenerate non-zero)')
    return dict(point=float(point), bootstrap_mean=mean, bias=bias, SEboot=se, b=b,
                bias_flag=bool(isinstance(b, float) and abs(b)>0.1))


def describe(values, ids):
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        return dict(estimable=False, mean='Infinity' if np.isposinf(values).any() else None,
                    reason='non-finite seed-window values; no ordinary SE/bootstrap/paired cost inference',
                    n=int(len(values)), finite_n=int(np.isfinite(values).sum()))
    boot = values[ids].mean(axis=1)
    point = float(np.mean(values))
    se = float(np.std(values, ddof=1)/np.sqrt(len(values)))
    return dict(estimable=True, n=len(values), mean=point, se=se,
                ci95=np.quantile(boot,[.025,.975],method='linear').tolist(),
                ci90=np.quantile(boot,[.05,.95],method='linear').tolist(),
                compatibility_ci=np.quantile(boot,[.05/12,1-.05/12],method='linear').tolist(),
                diagnostic=bias_diagnostic(point,boot))


def paired(a, b, ids):
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        return dict(estimable=False, reason='non-finite side: ordinary paired cost difference, SE and bootstrap prohibited')
    out = describe(np.asarray(a)-b, ids)
    out['delta_over_se'] = (out['mean']/out['se'] if out['se'] else
                            ('zero difference' if out['mean']==0 else 'nonzero difference / zero SE'))
    return out


def ratio(a, b, ids):
    if not np.isfinite(a).all() or not np.isfinite(b).all() or np.any(np.asarray(b)<=0):
        return dict(estimable=False, reason='nonpositive or non-finite denominator / non-finite numerator')
    with np.errstate(over='ignore',divide='ignore',invalid='ignore'):
        point = float(np.mean(a)/np.mean(b))
        boot = np.mean(np.asarray(a)[ids],axis=1)/np.mean(np.asarray(b)[ids],axis=1)
    if not np.isfinite(point) or not np.isfinite(boot).all():
        return dict(estimable=False, reason='non-finite ratio point or bootstrap; no draws removed')
    return dict(estimable=True, mean=point, ci90=np.quantile(boot,[.05,.95],method='linear').tolist(),
                ci95=np.quantile(boot,[.025,.975],method='linear').tolist(),diagnostic=bias_diagnostic(point,boot))


def compatibility(point, interval, target, tolerance):
    if interval is None or not np.isfinite(interval).all():
        return '不可估计'
    lo, hi = interval
    if hi < target-tolerance or lo > target+tolerance:
        return '不相容'
    return '相容' if abs(point-target)<=tolerance else '抽样不确定下仍相容，需披露'


def raw_m(ci):
    return bool(ci[0] > .15)


def raw_e(ci):
    return bool(ci[0] > -.02 and ci[1] < .02)


def raw_j(dw, dv):
    return bool(dw['estimable'] and dv['estimable'] and dw['ci90'][1]<0 and dv['ci90'][1]<0)


def gate(eligible, m, e, j, j_estimable=True):
    if not eligible:
        return {k:'未检验（质量门关闭）' for k in ('M','E','J')}
    states={'M':'通过' if m else '未通过'}
    states['E'] = ('通过' if e else '未通过') if m else '未检验（前序未通过）'
    states['J'] = ('通过' if j else '未通过') if m and e else '未检验（前序未通过）'
    if m and e and not j_estimable:
        states['J']='不可估计'
    return states


def conclusion(states, sensitive_m, sensitive_e_numeric, sensitive_e_algorithm, l2_complete):
    if states['M'].startswith('未检验'):
        m = Q+'基准/负对照统计相容性未通过，主推断被阻断；原始 M/E/J 只作描述。'
    elif states['M']=='未通过':
        m = Q+'预注册的超过 15 个百分点收益门槛未通过。'
    elif sensitive_m:
        m = Q+'主设定下预注册门槛通过，但该结论对积分/索引约定敏感。'
    elif not l2_complete:
        m = Q+'主设定下预注册收益门槛通过，预定稳健性检查未完成。'
    else:
        m = Q+'支持平均覆盖改善超过 15 个百分点。'
    e = None
    if states['M']=='通过':
        if states['E']=='未通过':
            e=Q+'覆盖等效未获支持。'
        elif states['E']=='通过':
            restrictions=[]
            if sensitive_e_algorithm:
                restrictions.append('等效结论依赖本稿预先固定的 A')
            if sensitive_e_numeric:
                restrictions.append('等效结论对积分/索引约定敏感')
            if not l2_complete:
                restrictions.append('约定敏感性检查未完成')
            e=Q+('主设定下等效门槛通过，但'+'，且'.join(restrictions)+'。' if restrictions else
                 '主设定下平均覆盖在 ±2 个百分点内等效。')
    return dict(M=m,E=e)
