#!/usr/bin/env python3
"""Standalone verifier: stdlib + NumPy only; never imports the main package.

Reconstructs every interval's scores, all seed metrics, bootstrap summaries,
prequential widths and decisions from saved raw observations and interval logs.
Writes a distinct audit artifact; cannot overwrite the sealed experiment.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[name]='1'
import argparse
import hashlib
import json
import time
from pathlib import Path
from statistics import NormalDist
import numpy as np

WINDOWS={'all':(0,800),'pre':(0,160),'early50':(160,210),'transition20':(160,180),'pure30':(180,210),'high':(160,800)}
CONDITIONS=['constant','shift']
METRICS=['coverage','width','winkler']
Z=NormalDist().inv_cdf(.975)


def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def close(a,b,metric='coverage'):
    if metric=='winkler': np.testing.assert_allclose(a,b,atol=1e-10,rtol=1e-10)
    else: np.testing.assert_allclose(a,b,atol=1e-12,rtol=0)


def bootstrap(x,ids):
    x=np.asarray(x)
    if not np.isfinite(x).all(): return None
    samples=np.mean(x[ids],axis=1)
    return dict(mean=float(x.mean()),se=float(x.std(ddof=1)/np.sqrt(len(x))),
                ci95=np.quantile(samples,[.025,.975],method='linear'),
                ci90=np.quantile(samples,[.05,.95],method='linear'),
                compatibility_ci=np.quantile(samples,[.05/12,1-.05/12],method='linear'),
                bootstrap_mean=float(samples.mean()),bias=float(samples.mean()-x.mean()),SEboot=float(samples.std(ddof=1)))


def compare_summary(x,expected,ids,metric):
    d=bootstrap(x,ids)
    if d is None:
        assert expected['estimable'] is False
        assert 'ci95' not in expected and 'se' not in expected
        return
    assert expected['estimable'] is True
    for field in ('mean','se','ci95','ci90','compatibility_ci'): close(d[field],expected[field],metric)
    diag=expected['diagnostic']
    for field in ('bootstrap_mean','bias','SEboot'): close(d[field],diag[field],metric)
    b=d['bias']/d['SEboot'] if d['SEboot'] else ('N/A (degenerate)' if d['bias']==0 else 'ERROR (degenerate non-zero)')
    if isinstance(b,str): assert diag['b']==b
    else: close(b,diag['b'],metric)
    assert diag['bias_flag']==bool(isinstance(b,float) and abs(b)>.1)


def features(lags):
    return np.column_stack([np.ones(len(lags)),lags]+[lags[:,a]*lags[:,b] for a in range(20) for b in range(a,20)])


def main():
    p=argparse.ArgumentParser();p.add_argument('directory',nargs='?',default='artifacts/v1');p.add_argument('--output',default='artifacts/recompute_v1.json')
    args=p.parse_args();root=Path(args.directory);dest=Path(args.output)
    if dest.exists(): raise RuntimeError('refuse to overwrite independent audit')
    start=time.monotonic();manifest=json.loads((root/'artifact_manifest.json').read_text())
    for name,item in manifest.items():
        assert digest(root/name)==item['sha256'],name
        assert (root/name).stat().st_size==item['bytes']
    result=json.loads((root/'result.json').read_text());record=json.loads((root/'run_record.json').read_text())
    assert record['status']=='complete'
    attempts=json.loads((root/'attempts.json').read_text())
    assert len(attempts)==300 and all(x['status']=='complete' and x['exit_code']==0 for x in attempts)
    ids=np.load(root/'bootstrap_ids.npy')
    np.testing.assert_array_equal(ids,np.random.Generator(np.random.PCG64(2026110899)).integers(0,100,(10000,100)))
    innovations=np.load(root/'innovations.npy')
    for r in range(100): np.testing.assert_array_equal(innovations[r],np.random.Generator(np.random.PCG64(2026100800+r)).standard_normal(3221))
    latent=np.load(root/'trajectories.npz')
    # Re-integrate in separate scalar code, with synchronous right-hand sides.
    for setup,h,stride in [('main',.005,20),('half_step',.0025,40)]:
        x=u=v=0.;track=[]
        for k in range(3721*stride):
            dx=np.sin(k*h)-.1*x+.5*np.tanh(u);du=x-.2*u-2*np.pi*v;dv=2*np.pi*u-.2*v
            x,u,v=x+h*dx,u+h*du,v+h*dv
            if (k+1)%stride==0 and (k+1)//stride>500: track.append((x,u,v))
        close(np.array(track),latent[setup])
    saved=np.load(root/'seed_metrics.npz');arrays={};interval_count=0
    adaptive={};fit_count=0
    for setup in ('main','half_step','index20'):
        names=result['methods'][setup];off=20 if setup=='index20' else 21
        target=np.arange(3200)+off;lags=target[:,None]-np.arange(1,21)
        saved_indices=np.load(root/setup/'indices.npz')
        np.testing.assert_array_equal(saved_indices['inputs'],lags);np.testing.assert_array_equal(saved_indices['targets'],target)
        a=np.empty((2,100,len(names),6,3));arrays[setup]=a
        for r in range(100):
            f=np.load(root/setup/f'{r:03d}_fit.npz');ys=np.load(root/setup/f'{r:03d}_observations.npz')
            scalar=latent['half_step' if setup=='half_step' else 'main'][:,0]
            design=features(ys['constant'][lags[:1600]])
            gram=design.T@design+16*np.eye(231);rhs=design.T@ys['constant'][target[:1600]]
            np.testing.assert_allclose(gram@f['w'],rhs,rtol=1e-10,atol=1e-9)
            np.testing.assert_allclose(f['chol']@f['chol'].T,gram,rtol=1e-12,atol=1e-10)
            tau=np.sqrt(np.mean((ys['constant'][target[:1600]]-design@f['w'])**2));close(tau,f['tau'])
            cal=ys['constant'][target[1600:2400]]-features(ys['constant'][lags[1600:2400]])@f['w']
            close(cal,f['cal_resid']);fit_count+=1
            df=231-16*np.sum(np.linalg.solve(f['chol'],np.eye(231))**2);close(df,f['df'])
            close(tau*np.sqrt(1600/(1600-df)),f['tau_df'])
            for ci,c in enumerate(CONDITIONS):
                sigma=np.full(3221,.15)
                if c=='shift': sigma[2560+off:3200+off]=.35
                np.testing.assert_array_equal(ys[c],scalar+sigma*innovations[r])
                z=np.load(root/setup/f'{r:03d}_{c}_intervals.npz')
                y=ys[c][target[2400:]];np.testing.assert_array_equal(y,z['target'])
                test_design=features(ys[c][lags[2400:]])
                np.testing.assert_allclose(test_design@f['w'],z['center'],atol=1e-12,rtol=0)
                solution=np.linalg.solve(f['chol'],test_design.T)
                lev=np.sum(solution**2,axis=0);close(lev,z['leverage'])
                low,up=z['lower'],z['upper']
                hits=(low<=y)&(y<=up);width=up-low
                wink=width.copy();finite=np.isfinite(width)
                for mi in range(len(names)):
                    ok=finite[mi]
                    wink[mi,ok]+=40*np.maximum(low[mi,ok]-y[ok],0)+40*np.maximum(y[ok]-up[mi,ok],0)
                np.testing.assert_array_equal(hits,z['hit']);close(width,z['width']);close(wink,z['winkler'],'winkler')
                close(low,z['center']-z['halfwidth']);close(up,z['center']+z['halfwidth'])
                residual=np.concatenate([cal,y-z['center']])
                alphas={m:.05 for m in names if m.startswith('A')}
                frozen_q=np.sort(np.abs(cal))[760]
                for t in range(800):
                    hist=residual[700+t:800+t]
                    rms=Z*np.sqrt(np.mean(hist**2))
                    expected={'B':Z*tau*np.sqrt(1+lev[t]),'F':frozen_q,'G':rms,
                              'G_L':rms*np.sqrt(1+lev[t]),'B_df':Z*float(f['tau_df'])*np.sqrt(1+lev[t]),
                              'rolling_F':np.sort(np.abs(hist))[95]}
                    for m in alphas:
                        mi=names.index(m);v=np.abs(residual[:800+t]) if m=='A_full' else np.abs(hist)
                        alpha=alphas[m]
                        close(alpha,z['alpha'][mi,t]);close(v.max(),z['history_max'][mi,t])
                        expected_first=-800 if m=='A_full' else t-100
                        assert z['history_first'][mi,t]==expected_first and z['history_last'][mi,t]==t-1
                        if m=='A_order':
                            k=int(np.ceil((len(v)+1)*(1-alpha)))
                            h=0 if k==0 else (np.inf if k==len(v)+1 else np.sort(v)[k-1])
                        elif m=='A_inf' and alpha==0: h=np.inf
                        else: h=np.quantile(v,1-alpha,method='linear')
                        expected[m]=h
                        lo,hi=(.001,.999) if m=='A_clip' else (0,1)
                        alphas[m]=float(np.clip(alpha+.01*(.05-int(not hits[mi,t])),lo,hi))
                    # Vectorized comparison avoids a Python assertion per method.
                    close([expected[m] for m in names],z['halfwidth'][:,t])
                    gi=names.index('G')
                    assert z['history_first'][gi,t]==t-100 and z['history_last'][gi,t]==t-1
                if 'G_L' in names:
                    assert np.all(hits[names.index('G_L')]>=hits[names.index('G')])
                for wi,(w,(lo,hi)) in enumerate(WINDOWS.items()):
                    for mi,m in enumerate(names):
                        a[ci,r,mi,wi]=[hits[mi,lo:hi].mean(),width[mi,lo:hi].mean(),wink[mi,lo:hi].mean()]
                        if m.startswith('A'):
                            key=f'{setup}/{c}/{m}/{w}'
                            d=adaptive.setdefault(key,dict(infinite_points=0,infinite_seeds=0,alpha_zero=0,alpha_lt_001=0,halfwidth_equals_history_max=0,point_denominator=100*(hi-lo),seed_denominator=100))
                            inf=np.isinf(width[mi,lo:hi]);np.testing.assert_array_equal(inf,z['infinite'][mi,lo:hi])
                            d['infinite_points']+=int(inf.sum());d['infinite_seeds']+=int(inf.any())
                            d['alpha_zero']+=int(np.sum(z['alpha'][mi,lo:hi]==0));d['alpha_lt_001']+=int(np.sum(z['alpha'][mi,lo:hi]<.01))
                            d['halfwidth_equals_history_max']+=int(np.sum(z['halfwidth'][mi,lo:hi]==z['history_max'][mi,lo:hi]))
                interval_count+=low.size
            if (r+1)%25==0: print(f'independent {setup}: {r+1}/100 ({time.monotonic()-start:.1f}s)',flush=True)
        for ki,metric in enumerate(METRICS): close(a[...,ki],saved[setup][...,ki],metric)
    assert adaptive==result['adaptive_diagnostics']
    def get(setup,c,m,w,metric): return arrays[setup][CONDITIONS.index(c),:,result['methods'][setup].index(m),list(WINDOWS).index(w),METRICS.index(metric)]
    for key,stats in result['summary'].items():
        setup,c,m,w=key.split('/')
        for metric in METRICS: compare_summary(get(setup,c,m,w,metric),stats[metric],ids,metric)
    for key,stats in result['comparisons'].items():
        name,c,w=key.split('/')
        if name.startswith(('half_step-','index20-')):
            setup,rest=name.split('-',1)
            if rest in ('G-B','G-A'): m,b=rest.split('-');base=setup
            else: m=b=rest;base='main'
        else: setup=base='main';m,b=name.split('-')
        for metric in METRICS:
            x=get(setup,c,m,w,metric);y=get(base,c,b,w,metric)
            if not np.isfinite(x).all() or not np.isfinite(y).all():
                assert stats[metric]['estimable'] is False and 'mean' not in stats[metric]
            else: compare_summary(x-y,stats[metric],ids,metric)
    for metric,j in result['J_statistics'].items():
        g=get('main','shift','G','high',metric);a=get('main','shift','A','high',metric)
        compare_summary(g-1.1*a,j['D'],ids,metric)
        ratios=g[ids].mean(1)/a[ids].mean(1);point=g.mean()/a.mean()
        close(j['ratio']['mean'],point,metric);close(j['ratio']['ci90'],np.quantile(ratios,[.05,.95],method='linear'),metric)
        close(j['ratio']['diagnostic']['bias'],ratios.mean()-point,metric)
    quality=True
    for item in result['N_F']['items']:
        d=bootstrap(get('main',item['condition'],item['method'],item['window'],'coverage'),ids)
        lo,hi=d['compatibility_ci'];c=item['target'];tol=item['tolerance']
        compatible=not (hi<c-tol or lo>c+tol);quality &= compatible
        state=('相容' if abs(d['mean']-c)<=tol else '抽样不确定下仍相容，需披露') if compatible else '不相容'
        assert item['state']==state;close(item['interval'],d['compatibility_ci'])
    m=bootstrap(get('main','shift','G','high','coverage')-get('main','shift','B','high','coverage'),ids)['ci95'][0]>.15
    eci=bootstrap(get('main','shift','G','high','coverage')-get('main','shift','A','high','coverage'),ids)['ci90']
    e=eci[0]>-.02 and eci[1]<.02
    j=all(result['J_statistics'][metric]['D']['ci90'][1]<0 for metric in ('width','winkler'))
    raw={'M':bool(m),'E':bool(e),'J':bool(j)}
    for k,v in raw.items(): assert result['raw'][k]==v
    if not quality: states={k:'未检验（质量门关闭）' for k in raw}
    else:
        states={'M':'通过' if m else '未通过','E':('通过' if e else '未通过') if m else '未检验（前序未通过）',
                'J':('通过' if j else '未通过') if m and e else '未检验（前序未通过）'}
    assert states==result['gated'] and quality==result['N_F']['passed']
    payload={'status':'passed','raw_artifact_manifest_sha256':digest(root/'artifact_manifest.json'),
             'result_sha256':digest(root/'result.json'),'intervals_recomputed':interval_count,'fits_checked':fit_count,
             'seed_attempts':len(attempts),'all_summary_and_comparison_statistics_checked':True,
             'all_adaptive_diagnostics_checked':True,'all_prequential_widths_checked':True,
             'independent_raw':raw,'independent_gated':states,'quality_gate':quality,
             'coverage_width_atol':1e-12,'winkler_atol':1e-10,'winkler_rtol':1e-10,
             'wall_seconds':time.monotonic()-start}
    dest.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(payload,ensure_ascii=False),flush=True)

if __name__=='__main__': main()
