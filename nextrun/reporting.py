"""Reporting supplements shared by tests and the report renderer."""
def main_paragraph(result):
    early=result['comparisons']['G-B/shift/early50']['coverage']
    text=result['conclusion']['M']
    text+=(f" 前 50 点 G−B={100*early['mean']:.3f} 个百分点，SE={100*early['se']:.3f}，"
           f"95% 配对区间 [{100*early['ci95'][0]:.3f}, {100*early['ci95'][1]:.3f}] 个百分点（无确认性资格）。")
    deviations=[]
    for m in ('B','G'):
        d=result['width_reference'][m]
        if d['over_10_percent']: deviations.append(f"{m} 相对差 {100*d['relative_difference']:+.2f}%")
    if deviations: text+=' 覆盖结论不附带宽度复现（'+'；'.join(deviations)+'）。'
    flags=[];unavailable=[]
    for metric,item in result['J_statistics'].items():
        for stat,value in item.items():
            d=value.get('diagnostic')
            if not d: unavailable.append(f'{metric}/{stat}: '+value.get('reason','不可估计'))
            elif d['bias_flag']: flags.append(f"{metric}/{stat}: b={d['b']:.4f}")
            elif isinstance(d['b'],str): unavailable.append(f"{metric}/{stat}: {d['b']}")
    text+=' J 门控状态：'+result['gated']['J']+'。'
    if flags: text+=' bootstrap 偏差相对 SE 超过预定阈值（'+'；'.join(flags)+'）；描述比值不参与线性 J 的确认判定。'
    else: text+=' J 两项线性差及描述比值均未触发 |b|>0.1 旗标；这不证明区间覆盖准确。'
    if unavailable: text+=' 诊断限制：'+'；'.join(unavailable)+'。'
    return text
