#!/usr/bin/env python3
"""Post-analysis integrity and presentation checks; no new experiments or thresholds."""
import os
for variable in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[variable]='1'
import numpy as np
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    r=json.loads((ROOT/'evaluation/result.json').read_text())
    audit=json.loads((ROOT/'evaluation/independent_recompute.json').read_text())
    assert audit['status']=='passed'
    assert sha(ROOT/'evaluation/result.json')==audit['result_sha256']
    assert sha(ROOT/'evaluation/artifact_manifest.json')==audit['raw_artifact_manifest_sha256']
    assert sha(ROOT/'PREREG.md')=='9d486db80eb5faca0d5aeb21c60ae80088991514f2360fd99208b22e71c77991'
    assert not subprocess.check_output(['git','diff','ed8d08bdd29e6084ccb42345896352ea1e62ad8e','--','PREREG.md'],cwd=ROOT)
    for line in (ROOT/'EXECUTION_SHA256SUMS').read_text().splitlines():
        digest,name=line.split('  ',1);assert sha(ROOT/name)==digest,name
    for name in ('result.json','run_record.json','attempts.json','fit_diagnostics.json','seed_metrics.npz'):
        assert sha(ROOT/'evaluation'/name)==sha(ROOT/'artifacts/v1'/name),name
    posts=re.findall(r'<!-- Post (\d+): (\d+) characters -->\n([^\n]+)',(ROOT/'x_post.md').read_text())
    assert 1<=len(posts)<=4
    for i,count,text in posts: assert int(count)==len(text)<=270
    report=(ROOT/'REPORT.md').read_text()
    assert '草稿：独立复算尚未完成' not in report
    assert r['conclusion']['M'] in report
    if not r['N_F']['passed']:
        assert all(s=='未检验（质量门关闭）' for s in r['gated'].values())
        assert '主推断被阻断' in report
    for phrase in ('无确认性资格','覆盖结论不附带宽度复现','bootstrap','首次公开精确时间待核实'):
        assert phrase in report
    for image in re.findall(r'!\[[^\]]*\]\(([^)]+)\)',report):assert (ROOT/image).is_file()
    for path in ('REPORT.md','README.md','EXECUTION.md'):
        for target in re.findall(r'\]\(([^)]+)\)',(ROOT/path).read_text()):
            if not target.startswith(('http:','https:','#')):
                assert (ROOT/target.split('#')[0]).exists(),(path,target)
    # Independently map numeric summaries to N, C95, exploratory gates and L2 flags.
    core={m:r['summary'][f'main/shift/{m}/high']['coverage'] for m in r['methods']['main']}
    expected_n=all(abs(core[m]['mean']-v)<=.02 for m,v in [('B',.726),('G',.929),('A',.938)])
    assert expected_n==r['N']['passed']
    for m,d in core.items():assert (d['ci90'][0]>.93 and d['ci90'][1]<.97)==r['C95'][m]
    for k,name in [('M_L','G_L-B'),('M_against_B_df','G-B_df')]:
        assert (r['comparisons'][f'{name}/shift/high']['coverage']['ci95'][0]>.15)==r['raw'][k]
    sm=se_num=se_alg=False
    for setup in ('half_step','index20'):
        cov={m:r['comparisons'][f'{setup}-{m}/shift/high']['coverage']['mean'] for m in ('B','G','A')}
        m_ci=r['comparisons'][f'{setup}-G-B/shift/high']['coverage']['ci95']
        e_ci=r['comparisons'][f'{setup}-G-A/shift/high']['coverage']['ci90']
        m=m_ci[0]>.15;e=e_ci[0]>-.02 and e_ci[1]<.02
        d=r['sensitivity']['details'][setup]
        assert d['coverage_changes']==cov and d['raw_M']==m and d['raw_E']==e
        sm |= any(abs(v)>.01 for v in cov.values()) or m!=r['raw']['M']
        se_num |= e!=r['raw']['E']
    seed_metrics=np.load(ROOT/'evaluation/seed_metrics.npz')['main']
    ids=np.load(ROOT/'artifacts/v1/bootstrap_ids.npy')
    methods=r['methods']['main'];high=list(r['windows']).index('high')
    for variant in ('A_full','A_clip','A_order','A_inf'):
        delta=seed_metrics[1,:,methods.index('G'),high,0]-seed_metrics[1,:,methods.index(variant),high,0]
        resampled=delta[ids].mean(axis=1)
        d=r['sensitivity']['details'][variant];lo,hi=d['G_minus_variant']['ci90']
        np.testing.assert_allclose([lo,hi],np.quantile(resampled,[.05,.95],method='linear'),atol=1e-12,rtol=0)
        e=lo>-.02 and hi<.02;assert e==d['raw_E']
        se_alg |= e!=r['raw']['E']
    assert [sm,se_num,se_alg]==[r['sensitivity'][k] for k in ('S_M','S_E_numeric','S_E_algorithm')]
    candidates=[p for p in ROOT.rglob('*') if p.is_file() and not any(x in p.relative_to(ROOT).parts for x in ('.git','.venv','artifacts','__pycache__'))]
    assert all(p.stat().st_size<2_000_000 for p in candidates),'file exceeds delivery size cap'
    result={'status':'passed','prereg_unchanged':True,'execution_freeze_checksums_passed':True,
            'sealed_result_matches_delivery':True,'report_links_and_figures_present':True,
            'main_and_sensitivity_decisions_verified':True,'post_characters':[len(t) for _,_,t in posts],
            'largest_delivery_file_bytes':max(p.stat().st_size for p in candidates),
            'branch':subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip()}
    assert result['branch']=='run/v1'
    (ROOT/'evaluation/delivery_validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))

if __name__=='__main__':main()
