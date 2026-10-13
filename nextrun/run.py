"""Sealed, serial execution. Refuses to overwrite an earlier analysis directory."""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
import json
import importlib.metadata
import sys
import resource
import subprocess
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from .core import trajectory, indices, noise_scales, observe, fit, feature, evaluate, A_METHODS
from .io import sha, write_json
from .analyze import analyze


def utc(): return datetime.now(timezone.utc).isoformat()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='artifacts/v1')
    args=parser.parse_args(); root=Path(args.output)
    if root.exists(): raise RuntimeError('refuse overwrite: use a new version and disclose any rerun')
    config=json.loads(Path('config.json').read_text())
    for line in Path('EXECUTION_SHA256SUMS').read_text().splitlines():
        digest,name=line.split('  ',1)
        assert sha(name)==digest, f'freeze mismatch: {name}'
    assert sha('PREREG.md')==config['protocol_sha256']
    software=json.loads(Path('software_lock.json').read_text())
    assert sys.version.split()[0]==software['python']
    for package,version in software['packages'].items():
        assert importlib.metadata.version(package)==version, f'dependency changed: {package}'
    dirty=subprocess.check_output(['git','status','--porcelain'],text=True)
    assert not dirty, 'execution must start from a clean committed freeze'
    head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    start=time.monotonic();root.mkdir(parents=True)
    resource.setrlimit(resource.RLIMIT_AS,(config['memory_bytes'],config['memory_bytes']))
    resource.setrlimit(resource.RLIMIT_CPU,(7100,7150))
    # SIGTERM/timeout and Python failures leave partial logs; never upgrade them to complete.
    record={'started_utc':utc(),'execution_commit':head,'protocol_commit':config['protocol_commit'],
            'protocol_sha256':config['protocol_sha256'],'status':'running','configuration':config}
    write_json(root/'run_record.json',record)
    attempts=[dict(setup=s,replicate=r,seed=config['seed_start']+r,status='not_started',exit_code=None,stage=None,seconds=None)
              for s in config['setups'] for r in range(100)]
    write_json(root/'attempts.json',attempts)
    def budget():
        if time.monotonic()-start>7000: raise TimeoutError('two-hour budget stop, allowing finalization margin')
        if sum(f.stat().st_size for f in root.rglob('*') if f.is_file())>config['artifact_bytes']:
            raise RuntimeError('artifact disk budget exceeded')
    try:
        ids=np.random.Generator(np.random.PCG64(config['bootstrap_seed'])).integers(0,100,size=(10000,100))
        np.save(root/'bootstrap_ids.npy',ids)
        innovations=np.array([np.random.Generator(np.random.PCG64(config['seed_start']+r)).standard_normal(3221) for r in range(100)])
        np.save(root/'innovations.npy',innovations)
        main_path=trajectory(config['h']);half_path=trajectory(config['half_h'])
        np.savez_compressed(root/'trajectories.npz',main=main_path,half_step=half_path,
                            physical_time=(501+np.arange(3221))*.1)
        trajectory_difference={'max_abs_x':float(np.max(np.abs(main_path[:,0]-half_path[:,0]))),
                               'rms_x':float(np.sqrt(np.mean((main_path[:,0]-half_path[:,0])**2))),
                               'max_abs_all_states':float(np.max(np.abs(main_path-half_path)))}
        arrays={};methods={};diagnostics={};fit_diagnostics=[];curves={}
        for setup in config['setups']:
            methods[setup]=config['methods'] if setup=='main' else config['methods'][:4]
            method_list=methods[setup]; directory=root/setup;directory.mkdir()
            offset=20 if setup=='index20' else 21
            path=half_path if setup=='half_step' else main_path
            inputs, targets=indices(offset)
            np.savez_compressed(directory/'indices.npz',inputs=inputs,targets=targets,train=np.arange(1600),
                                calibration=np.arange(1600,2400),test=np.arange(2400,3200),test_index=np.arange(800))
            arr=np.empty((2,100,len(method_list),6,3)); arrays[setup]=arr
            curves[setup]=np.zeros((2,len(method_list),800,3))
            diag={}
            for r in range(100):
                budget(); a=next(x for x in attempts if x['setup']==setup and x['replicate']==r)
                seed_start=time.monotonic();a.update(status='running',stage='observations')
                write_json(root/'attempts.json',attempts)
                ys={c:observe(path,innovations[r],noise_scales(offset,c)) for c in config['conditions']}
                assert np.array_equal(ys['constant'][:2400+offset],ys['shift'][:2400+offset])
                np.savez_compressed(directory/f'{r:03d}_observations.npz',**ys)
                a['stage']='fit'
                model=fit(ys['constant'][inputs[:1600]],ys['constant'][targets[:1600]])
                cal_resid=ys['constant'][targets[1600:2400]]-feature(ys['constant'][inputs[1600:2400]])@model.w
                np.savez_compressed(directory/f'{r:03d}_fit.npz',w=model.w,chol=model.chol,cal_resid=cal_resid,
                                    tau=model.tau,df=model.df,tau_df=model.tau_df)
                fit_diagnostics.append(dict(setup=setup,replicate=r,df=model.df,tau=model.tau,
                                            tau_df=model.tau_df,scale_ratio=model.tau_df/model.tau))
                for ci,condition in enumerate(config['conditions']):
                    a['stage']='predict_'+condition
                    log=evaluate(model,ys[condition],targets[2400:],cal_resid,method_list)
                    np.savez_compressed(directory/f'{r:03d}_{condition}_intervals.npz',methods=np.array(method_list),**log)
                    for wi,(window,(lo,hi)) in enumerate(config['windows'].items()):
                        for mi,method in enumerate(method_list):
                            for ki,metric in enumerate(('hit','width','winkler')):
                                arr[ci,r,mi,wi,ki]=np.mean(log[metric][mi,lo:hi])
                            if method in A_METHODS:
                                key=f'{setup}/{condition}/{method}/{window}'
                                d=diag.setdefault(key,dict(infinite_points=0,infinite_seeds=0,alpha_zero=0,alpha_lt_001=0,
                                                           halfwidth_equals_history_max=0,point_denominator=100*(hi-lo),seed_denominator=100))
                                inf=log['infinite'][mi,lo:hi]
                                d['infinite_points']+=int(inf.sum());d['infinite_seeds']+=int(inf.any())
                                d['alpha_zero']+=int(np.sum(log['alpha'][mi,lo:hi]==0))
                                d['alpha_lt_001']+=int(np.sum(log['alpha'][mi,lo:hi]<.01))
                                d['halfwidth_equals_history_max']+=int(np.sum(log['halfwidth'][mi,lo:hi]==log['history_max'][mi,lo:hi]))
                    for ki,metric in enumerate(('hit','width','winkler')):
                        curves[setup][ci,:,:,ki]+=log[metric]/100
                a.update(status='complete',exit_code=0,stage='complete',seconds=time.monotonic()-seed_start)
                write_json(root/'attempts.json',attempts)
                if (r+1)%10==0: print(f'{setup}: {r+1}/100 complete ({time.monotonic()-start:.1f}s)',flush=True)
            diagnostics.update(diag)
        budget()
        np.savez_compressed(root/'seed_metrics.npz',**arrays)
        np.savez_compressed(root/'curves.npz',**curves)
        write_json(root/'fit_diagnostics.json',fit_diagnostics)
        result=analyze(arrays,methods,config,ids,diagnostics)
        result.update(trajectory_difference=trajectory_difference,execution_commit=head,
                      protocol_sha256=config['protocol_sha256'],methods=methods,windows=config['windows'])
        write_json(root/'result.json',result)
        record.update(status='complete',finished_utc=utc(),wall_seconds=time.monotonic()-start,
                      peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                      attempts_complete=sum(x['status']=='complete' for x in attempts))
        write_json(root/'run_record.json',record)
        manifest={str(p.relative_to(root)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(root.rglob('*')) if p.is_file()}
        write_json(root/'artifact_manifest.json',manifest)
        print('SEALED complete; results not interpreted by execution runner.',flush=True)
    except BaseException as exc:
        for a in attempts:
            if a['status']=='running': a.update(status='failed',exit_code=1,error=repr(exc))
        write_json(root/'attempts.json',attempts)
        record.update(status='failed',finished_utc=utc(),error=repr(exc),traceback=traceback.format_exc(),wall_seconds=time.monotonic()-start)
        write_json(root/'run_record.json',record)
        raise

if __name__=='__main__': main()
