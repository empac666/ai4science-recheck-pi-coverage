"""Aggregate all preregistered settings without selecting variants."""
import numpy as np
from .statistics import describe, paired, ratio, compatibility, raw_m, raw_e, raw_j, gate, conclusion

METRICS=('coverage','width','winkler')


def analyze(arrays, methods, config, ids, diagnostics):
    windows=list(config['windows'])
    out={'schema_version':1,'scope':'observation noise conditional on frozen x*, not a trajectory population',
         'summary':{},'comparisons':{},'adaptive_diagnostics':diagnostics}
    def values(setup,condition,method,window,metric):
        return arrays[setup][config['conditions'].index(condition),:,methods[setup].index(method),windows.index(window),METRICS.index(metric)]
    for setup in config['setups']:
        for condition in config['conditions']:
            for method in methods[setup]:
                for window in windows:
                    key='/'.join((setup,condition,method,window))
                    out['summary'][key]={metric:describe(values(setup,condition,method,window,metric),ids) for metric in METRICS}
    def compare(name,s1,m1,s0,m0):
        for c in config['conditions']:
            for w in windows:
                out['comparisons']['/'.join((name,c,w))]={metric:paired(values(s1,c,m1,w,metric),values(s0,c,m0,w,metric),ids) for metric in METRICS}
    for a,b in [('G','B'),('G','A'),('A_full','A'),('A_clip','A'),('A_order','A'),('A_inf','A'),
                ('G_L','B'),('rolling_F','F'),('B_df','B'),('G','B_df')]:
        compare(a+'-'+b,'main',a,'main',b)
    for setup in ('half_step','index20'):
        for method in ('B','F','G','A'):
            compare(setup+'-'+method,setup,method,'main',method)
        compare(setup+'-G-B',setup,'G',setup,'B')
        compare(setup+'-G-A',setup,'G',setup,'A')
    summary=out['summary']; comps=out['comparisons']
    core={m:summary['main/shift/'+m+'/high']['coverage'] for m in ('B','G','A')}
    targets={'B':.726,'G':.929,'A':.938}
    out['N']={'passed':all(abs(core[m]['mean']-targets[m])<=.02 for m in targets),
              'differences':{m:core[m]['mean']-targets[m] for m in targets}}
    nf=[]
    for c,w,m,target,tol in [('shift','high','B',.726,.05),('shift','high','F',.595,.02),
                            ('constant','all','B',.941,.02),('constant','all','F',.949,.02),
                            ('constant','all','G',.947,.02),('constant','all','A',.949,.02)]:
        d=summary[f'main/{c}/{m}/{w}']['coverage']
        nf.append(dict(condition=c,window=w,method=m,target=target,tolerance=tol,
                       point=d['mean'],difference=d['mean']-target,se=d['se'],
                       interval=d['compatibility_ci'], band=[target-tol,target+tol],
                       state=compatibility(d['mean'],d['compatibility_ci'],target,tol)))
    eligible=all(d['state'] not in ('不相容','不可估计') for d in nf)
    out['N_F']={'passed':eligible,'items':nf,'quantile_probabilities':[.05/12,1-.05/12]}
    j={}
    for metric in ('width','winkler'):
        a=values('main','shift','G','high',metric); b=values('main','shift','A','high',metric)
        ok=np.isfinite(a).all() and np.isfinite(b).all() and np.all(b>0)
        j[metric]={'D':describe(a-1.1*b,ids) if ok else {'estimable':False,'reason':'invalid J input'},
                   'ratio':ratio(a,b,ids)}
    mr=raw_m(comps['G-B/shift/high']['coverage']['ci95'])
    er=raw_e(comps['G-A/shift/high']['coverage']['ci90'])
    jr=raw_j(j['width']['D'],j['winkler']['D'])
    je=all(j[x]['D']['estimable'] for x in j)
    out['J_statistics']=j
    out['raw']={'M':mr,'E':er,'J':jr,'M_L':raw_m(comps['G_L-B/shift/high']['coverage']['ci95']),
                'M_against_B_df':raw_m(comps['G-B_df/shift/high']['coverage']['ci95'])}
    out['gated']=gate(eligible,mr,er,jr,je)
    out['eligibility']={'all_100_main_seeds':True,'implementation_checks_passed':True,'quality_gate':eligible}
    # This is deliberately AFTER the immutable main decisions above.
    sensitive_m=False; sensitive_e_num=False; sensitive_e_alg=False
    sensitivities={}
    for setup in ('half_step','index20'):
        deltas={m:comps[f'{setup}-{m}/shift/high']['coverage']['mean'] for m in ('B','G','A')}
        rm=raw_m(comps[f'{setup}-G-B/shift/high']['coverage']['ci95'])
        re=raw_e(comps[f'{setup}-G-A/shift/high']['coverage']['ci90'])
        sensitive_m |= any(abs(v)>.01 for v in deltas.values()) or rm != mr
        sensitive_e_num |= re != er
        sensitivities[setup]={'coverage_changes':deltas,'raw_M':rm,'raw_E':re,'M_flip':rm!=mr,'E_flip':re!=er}
    for variant in ('A_full','A_clip','A_order','A_inf'):
        d=paired(values('main','shift','G','high','coverage'),values('main','shift',variant,'high','coverage'),ids)
        re=raw_e(d['ci90']); sensitive_e_alg |= re!=er
        sensitivities[variant]={'G_minus_variant':d,'raw_E':re,'E_flip':re!=er}
    out['sensitivity']={'S_M':bool(sensitive_m),'S_E_numeric':bool(sensitive_e_num),
                        'S_E_algorithm':bool(sensitive_e_alg),'L2_complete':True,'details':sensitivities}
    out['conclusion']=conclusion(out['gated'],sensitive_m,sensitive_e_num,sensitive_e_alg,True)
    out['C95']={m:bool(summary[f'main/shift/{m}/high']['coverage']['ci90'][0]>.93 and
                      summary[f'main/shift/{m}/high']['coverage']['ci90'][1]<.97) for m in methods['main']}
    ref_width={'B':.858,'F':.645,'G':1.458,'A':1.544}
    out['width_reference']={m:{'reference':v,'relative_difference':summary[f'main/shift/{m}/high']['width']['mean']/v-1,
                             'over_10_percent':abs(summary[f'main/shift/{m}/high']['width']['mean']/v-1)>.1} for m,v in ref_width.items()}
    out['author_se_comparison']={m:{'author_rounded_se':se,'our_se':core[m]['se'],
                                   'difference_over_our_se':out['N']['differences'][m]/core[m]['se'],
                                   'our_se_over_author_se':core[m]['se']/se if se else None,
                                   'note':'author rounded zero is not zero uncertainty' if se==0 else 'scale comparison only'}
                                  for m,se in [('B',.001),('G',.001),('A',0.0)]}
    # Same IDs: a stronger nested interval cannot have a smaller bootstrap difference.
    g=values('main','shift','G','high','coverage'); gl=values('main','shift','G_L','high','coverage')
    b=values('main','shift','B','high','coverage')
    assert np.all(gl>=g) and np.all((gl-b)[ids].mean(1)>=(g-b)[ids].mean(1))
    assert comps['G_L-B/shift/high']['coverage']['ci95'][0]>=comps['G-B/shift/high']['coverage']['ci95'][0]
    return out
