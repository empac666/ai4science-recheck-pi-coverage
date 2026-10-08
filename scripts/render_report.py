#!/usr/bin/env python3
"""Presentation only: reads sealed numerical results, never changes decisions."""
import os
os.environ['MPLBACKEND']='Agg'
import csv
import json
import shutil
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from nextrun.reporting import main_paragraph
from nextrun.io import write_json, sha

ROOT=Path('artifacts/v1')
EVAL=Path('evaluation');EVAL.mkdir(exist_ok=True)
FIG=Path('figures');FIG.mkdir(exist_ok=True)
r=json.loads((ROOT/'result.json').read_text())
run=json.loads((ROOT/'run_record.json').read_text())
draft='--draft' in sys.argv
audit_path=Path('artifacts/recompute_v1.json')
audit=json.loads(audit_path.read_text()) if audit_path.exists() else {'status':'pending','intervals_recomputed':0,'wall_seconds':0}
assert run['status']=='complete' and (audit['status']=='passed' or draft)
assert not r['C95']['G'], 'review fixed narrative if rendering a different run'
for name in ('result.json','run_record.json','artifact_manifest.json','attempts.json','fit_diagnostics.json','seed_metrics.npz'):
    shutil.copyfile(ROOT/name,EVAL/name)
if audit_path.exists(): shutil.copyfile(audit_path,EVAL/'independent_recompute.json')


def f(x,d=6):
    if x is None:return 'N/A'
    if isinstance(x,str): return '∞' if x=='Infinity' else x
    return f'{x:.{d}f}'


def ci(d,key='ci95',scale=1):
    return '['+', '.join(f(scale*x) for x in d[key])+']' if d.get('estimable',True) else '不可估计'


def stat(d):
    return f(d.get('mean'))+(f" ± {f(d['se'])}" if 'se' in d else '')


def s(m,w='high',c='shift',setup='main'):return r['summary'][f'{setup}/{c}/{m}/{w}']

def cmp(name,w='high',c='shift'):return r['comparisons'][f'{name}/{c}/{w}']

# Complete machine-readable tabular companion, including every mandatory bias diagnostic.
with (EVAL/'bootstrap_diagnostics.csv').open('w',newline='') as stream:
    writer=csv.writer(stream,lineterminator=chr(10));writer.writerow(['family','key','metric','estimable','point','bootstrap_mean','bias','SEboot','b','bias_flag','reason'])
    for family in ('summary','comparisons'):
        for key,metrics in r[family].items():
            for metric,item in metrics.items():
                d=item.get('diagnostic',{})
                writer.writerow([family,key,metric,item['estimable'],d.get('point'),d.get('bootstrap_mean'),d.get('bias'),d.get('SEboot'),d.get('b'),d.get('bias_flag'),item.get('reason','')])
    for metric,items in r['J_statistics'].items():
        for name,item in items.items():
            d=item.get('diagnostic',{})
            writer.writerow(['J',name,metric,item['estimable'],d.get('point'),d.get('bootstrap_mean'),d.get('bias'),d.get('SEboot'),d.get('b'),d.get('bias_flag'),item.get('reason','')])

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':140})
colors={'B':'#4169a1','F':'#929292','G':'#d6761e','A':'#188466'}
curves=np.load(ROOT/'curves.npz')['main']
fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
for ci0,c in enumerate(('constant','shift')):
    ax=axes[ci0]
    for mi,m in enumerate(('B','F','G','A')):
        line=curves[ci0,mi,:,0]
        ax.plot(np.arange(800),line,alpha=.15,lw=.65,color=colors[m])
        # Partial windows before t=99 are deliberately not called 100-point estimates.
        tail=np.convolve(line,np.ones(100)/100,mode='valid')
        ax.plot(np.arange(99,800),tail,label=m,color=colors[m],lw=1.6)
    for edge in (160,180,210,260):ax.axvline(edge,color='black',ls=':',alpha=.35,lw=.8)
    ax.axhline(.95,color='black',ls='--',lw=.9)
    ax.set_ylim(.4,1.015);ax.set_ylabel('Coverage');ax.set_title(c+' noise: per-step mean (faint), trailing 100 (solid)',pad=27)
    boundaries=ax.secondary_xaxis('top');boundaries.set_xticks([160,180,210,260]);boundaries.tick_params(labelsize=8,labelrotation=45)
axes[0].legend(ncol=4);axes[1].set_xlabel('Zero-based test index')
fig.tight_layout();fig.savefig(FIG/'coverage_curves.png');plt.close(fig)

fig,axes=plt.subplots(2,3,figsize=(13,7))
for row,c in enumerate(('constant','shift')):
    w='all' if c=='constant' else 'high'
    for col,metric in enumerate(('coverage','width','winkler')):
        ax=axes[row,col]
        for k,m in enumerate(('B','F','G','A')):
            d=s(m,w,c)[metric];mean=d['mean'];lo,hi=d['ci95']
            ax.errorbar(k,mean,yerr=np.array([[mean-lo],[hi-mean]]),fmt='o',color=colors[m],capsize=4)
        ax.set_xticks(range(4),['B','F','G','A']);ax.set_title(c+': '+metric)
        if metric=='coverage':ax.axhline(.95,color='gray',ls='--',lw=.8)
        ax.set_ylabel('Mean and paired-ID bootstrap 95% CI')
fig.tight_layout();fig.savefig(FIG/'coverage_costs.png');plt.close(fig)

names=['half_step-B','half_step-G','half_step-A','index20-B','index20-G','index20-A','A_full-A','A_clip-A','A_order-A','A_inf-A','G_L-B','rolling_F-F','B_df-B']
fig,axes=plt.subplots(1,2,figsize=(12,7),sharey=True)
for ax,w in zip(axes,('early50','high')):
    for k,name in enumerate(names):
        d=cmp(name,w)['coverage'];mean=100*d['mean'];lo,hi=np.array(d['ci95'])*100
        ax.errorbar(mean,k,xerr=np.array([[mean-lo],[hi-mean]]),fmt='o',color='#4169a1',capsize=3)
    ax.axvline(0,color='gray',ls='--');ax.set_yticks(range(len(names)),names);ax.set_xlabel('Paired coverage change (percentage points)');ax.set_title(w+' (descriptive)')
axes[0].invert_yaxis();fig.tight_layout();fig.savefig(FIG/'sensitivity.png');plt.close(fig)

lines=['# NEXTrun：NGRC 预测区间独立复核报告','',
       '**执行与全部预定 L2 已完成；研究结论由锁定判据决定。** 本报告是对已知原论文表格差距的独立实现检验，不是盲分析，也不是作者代码或原始样本的逐位复现。',
       '',main_paragraph(r),'']
if r['conclusion']['E']:lines.extend([r['conclusion']['E'],''])
lines.extend(['核心数字'+('在预定 ±0.02 容差内接近原表。' if r['N']['passed'] else '未全部接近原表。')+
              ' 本报告的所有区间只度量冻结轨道 x* 上独立观测噪声重复的 Monte Carlo 不确定性；不包括初值、轨道或其他动力系统的泛化不确定性。',
              '', '## 版本、封存和公开记录','',
              f"- 协议锁定提交：`{run['protocol_commit']}`；PREREG SHA-256：`{run['protocol_sha256']}`，原件未改。",
              f"- 结果前执行冻结提交：`{run['execution_commit']}`；源文件清单见 [EXECUTION_SHA256SUMS](EXECUTION_SHA256SUMS)。",
              f"- 运行 UTC：`{run['started_utc']}` 至 `{run['finished_utc']}`；墙钟 {run['wall_seconds']:.2f} 秒，峰值 RSS {run['peak_rss_bytes']/2**20:.2f} MiB。300/300 次拟合尝试完成，单进程、BLAS 单线程，无 GPU。",
              f"- 封存结果 SHA-256：`{sha(ROOT/'result.json')}`；[完整结果 JSON](evaluation/result.json)、[全部产物哈希](evaluation/artifact_manifest.json)、[尝试清单](evaluation/attempts.json)。",
              '- [执行前公开读取记录](evidence/public_lock.json)确认匿名可读且协议哈希相同；GitHub events 接口 403，首次公开精确时间待核实。读取时间不是首次 push 时间；未使用本地 Git 时间冒充公开时间戳。',
              '- [偏离登记](DEVIATIONS.md)：用户指定的两小时硬限、首次公开时间证据限制，以及复用固定系统依赖的环境隔离限制；方法、种子、窗口与统计门槛保持预注册约定。',
              '', '## 实现、验收与独立复算','',
              '主实现采用同步 Euler、20 个带噪滞后、231 列特征和含截距惩罚的岭回归（Gram 加 16I）。两条件共享创新与拟合；测试期每次先固定所有区间，再揭示 y。预测和历史接口接收 owning copy 且只读，更新返回新状态。A 的主约定为 100 残差线性分位数、[0,1] 投影，不冒称为已核实的作者实现。',
              '', '11 项结果前验收通过，覆盖特征/惩罚、同步积分、观测映射与索引、历史隔离/防泄漏、端点/Winkler、A 更新、质量门、线性 J 和结论优先级。见 [验收记录](validation_before_run.txt) 与 [测试](tests/test_protocol.py)。',
              '',f"独立 [recompute.py](scripts/recompute.py) 不导入主实现，重积分两条轨道、复核 300 次拟合，重新检查 {audit['intervals_recomputed']:,} 个逐点区间的预测前半宽、全部覆盖/宽度/Winkler、种子级汇总、bootstrap 和门控；耗时 {audit['wall_seconds']:.2f} 秒。覆盖和宽度 atol=1e−12，Winkler atol=rtol=1e−10，全部通过。这里“独立”指独立代码路径与公式重写，不是独立团队审稿。见 [独立复算记录](evaluation/independent_recompute.json)。",
              '', '## 主结果与原表参照','',
              '下表为突变条件的高噪声窗口 {160,…,799}；每个种子先汇总再等权平均，SE 的统计单位为 100 个噪声重复。作者数字来自 [论文 v1 表 1](https://arxiv.org/html/2609.32169v1#S5.SS4)，不是本次结果；作者 A 的窗口、分位数与投影约定仍不唯一。',
              '', '| 方法 | 原表覆盖 | 本次覆盖 ± SE | 95% 区间 | 覆盖差（百分点） | 均宽 ± SE | Winkler ± SE | 宽度相对原表差 |',
              '|---|---:|---:|---|---:|---:|---:|---:|'])
for m,ref in [('B',.726),('F',.595),('G',.929),('A',.938)]:
    d=s(m);lines.append(f"| {m} | {ref} | {stat(d['coverage'])} | {ci(d['coverage'])} | {100*(d['coverage']['mean']-ref):+.3f} | {stat(d['width'])} | {stat(d['winkler'])} | {100*r['width_reference'][m]['relative_difference']:+.2f}% |")
lines.extend(['','原文没有 Winkler 参照值。原文 B/G 高噪声 SE 均舍入为 0.001，A 为 0.000；舍入零不意味着没有不确定性。',
              '', '| 方法 | 本次 SE / 作者显示 SE | 点差 / 本次 SE |', '|---|---:|---:|'])
for m,d in r['author_se_comparison'].items():lines.append(f"| {m} | {f(d['our_se_over_author_se'])} | {f(d['difference_over_our_se'])} |")
lines.extend(['','![覆盖和代价](figures/coverage_costs.png)','',
              '图 1：常噪声全测试段与突变条件高噪声全段的覆盖、均宽和 Winkler；误差线为按重复 ID 配对 bootstrap 的 95% 百分位区间。较高覆盖、接近 95% 与较低代价是不同判断。',
              '', '## 质量门与固定序列','',
              '六项 N-F 采用同一 10,000×100 重采样 ID 表，分位概率 0.05/12 与 1−0.05/12，端点相交仍相容。作者目标值仅作固定参照；这不是两篇研究的总体相等检验。',
              '', '| 条件/窗口/方法 | 点差 | SE | Bonferroni 相容区间 | 容差带 | 状态 |', '|---|---:|---:|---|---|---|'])
for d in r['N_F']['items']:
    lines.append(f"| {d['condition']}/{d['window']}/{d['method']} | {f(d['difference'])} | {f(d['se'])} | {d['interval']} | {d['band']} | {d['state']} |")
lines.extend(['',f"N：{'通过' if r['N']['passed'] else '未通过'}；N-F：{'通过' if r['N_F']['passed'] else '未通过'}。资格与主设置门控状态保持独立，不被 L2 结果回写。",
              '', '| 判据 | 点估计 | 预定区间 | 原始数值判据 | 门控状态 |', '|---|---:|---|---|---|'])
for label,name,key in [('M','G-B','ci95'),('E','G-A','ci90')]:
    d=cmp(name)['coverage'];lines.append(f"| {label} ({name}) | {f(d['mean'])} | {ci(d,key)} | {r['raw'][label]} | {r['gated'][label]} |")
for metric,item in r['J_statistics'].items():
    d=item['D'];lines.append(f"| J ({metric}: G−1.10A) | {f(d['mean'])} | {ci(d,'ci90')} | {r['raw']['J']}（联合） | {r['gated']['J']} |")
lines.extend(['','门槛使用未舍入值：M 的 95% 下界 >0.15；E 的 90% 区间完全在 (−0.02,0.02)；J 两项 90% 上界都 <0。质量门关闭时，原始布尔值不能支持主假设；失败不能解释为没有任何收益或论文错误。',
              '', '| J 量/描述比值 | 点估计 | 90% 区间 | bootstrap 均值 | bias | SEboot | b | 旗标 |', '|---|---:|---|---:|---:|---:|---:|---|'])
for metric,items in r['J_statistics'].items():
    for name,d in items.items():
        z=d['diagnostic'];lines.append(f"| {metric}/{name} | {f(d['mean'])} | {ci(d,'ci90')} | {f(z['bootstrap_mean'])} | {f(z['bias'],9)} | {f(z['SEboot'])} | {f(z['b'])} | {z['bias_flag']} |")
lines.extend(['','描述比值采用跨种子均值之比，不参与 J。所有方法/窗口 Winkler、预定有限 Winkler 差、D 和比值的诊断均见 [完整诊断 CSV](evaluation/bootstrap_diagnostics.csv)；非法无穷代价差未计算。',
              '', 'G 未满足探索性的接近名义覆盖判据。探索性 C95（高噪声 90% 区间严格在 (0.93,0.97) 内）：'+ '；'.join(m+'='+str(r['C95'][m]) for m in ('B','F','G','A'))+'。这不是证明真实覆盖至少 95%。',
              '', '## 时间窗口与适应过程','',
              '160 时目标首次高噪声而全部滞后仍低噪声；179 时 19 个高噪声滞后；180 起全部 20 个滞后高噪声。G/主 A 在 160/180/260 的残差 ID 分别为 60–159、80–179、160–259；高噪声主窗口的 540/640 个时刻有全高噪声残差历史。这不保证稳态，也不能用窗口比较因果分离两种过渡。',
              '', '![逐步平均覆盖及尾随覆盖](figures/coverage_curves.png)','',
              '图 2：淡线为 100 次重复逐测试步平均命中，实线为完整 100 点尾随平均，虚线标记 160/180/210/260 与名义覆盖。只描述适应过程，不逐步做显著性检验。',
              '', '| 条件 | 窗口 | 方法 | 覆盖 ± SE | 均宽 ± SE | Winkler ± SE |', '|---|---|---|---:|---:|---:|'])
for c in ('constant','shift'):
    for w in r['windows']:
        for m in ('B','F','G','A'):
            d=s(m,w,c);lines.append(f"| {c} | {w} | {m} | {stat(d['coverage'])} | {stat(d['width'])} | {stat(d['winkler'])} |")
lines.extend(['','## 全部预定 L2','',
              f"步长减半轨道 x 的最大绝对差 {r['trajectory_difference']['max_abs_x']:.9g}，RMS 差 {r['trajectory_difference']['rms_x']:.9g}；这里比较的是规定数值路径变体，不能假称同一 x*，也不建立轨道总体。",'',
              f"敏感旗标：S_M={r['sensitivity']['S_M']}，S_E数值={r['sensitivity']['S_E_numeric']}，S_E算法={r['sensitivity']['S_E_algorithm']}；全部 L2 完成。以下原始门槛仅描述，不重开确认序列。",'',
              '| 积分/索引变体 | B/G/A 覆盖变化 | 原始 M（是否翻转） | 原始 E（是否翻转） |','|---|---|---|---|'])
for setup in ('half_step','index20'):
    d=r['sensitivity']['details'][setup]
    delta='；'.join(f'{m}: {100*v:+.4f} pp' for m,v in d['coverage_changes'].items())
    lines.append(f"| {setup} | {delta} | {d['raw_M']} ({d['M_flip']}) | {d['raw_E']} ({d['E_flip']}) |")
lines.extend(['','| 比较（突变高噪声） | 覆盖配对差 ± SE | 95% 区间 | Δ/SE | 均宽差 ± SE | Winkler 差 ± SE |','|---|---:|---|---:|---:|---:|'])
for name in names+['G-B_df']:
    d=cmp(name);lines.append(f"| {name} | {stat(d['coverage'])} | {ci(d['coverage'])} | {f(d['coverage']['delta_over_se'])} | {stat(d['width']) if d['width']['estimable'] else '不可估计（含∞）'} | {stat(d['winkler']) if d['winkler']['estimable'] else '不可估计（含∞）'} |")
lines.extend(['','所有条件与窗口的种子内配对差、SE、95% 区间、Δ/SE 及代价偏差诊断见完整 JSON，不对短窗口套用 1 个百分点工程容差。',
              '', '![预定敏感性](figures/sensitivity.png)','',
              '图 3：预定比较在前 50 点及高噪声全段的配对覆盖变化和 95% 区间。不同变体含不同机制，不按接近原表程度选取替代主结果。',
              '', '| A 单因素变体 | G−变体覆盖差 | 90% 区间 | 原始 E | 相对主 E 翻转 |','|---|---:|---|---|---|'])
for m in ('A_full','A_clip','A_order','A_inf'):
    d=r['sensitivity']['details'][m];lines.append(f"| {m} | {f(d['G_minus_variant']['mean'])} | {ci(d['G_minus_variant'],'ci90')} | {d['raw_E']} | {d['E_flip']} |")
lines.extend(['','A_full 只改全部校准及已揭示历史；A_clip 只改投影 [0.001,0.999]；A_order 只改有限样本顺序统计量；A_inf 只改 α=0 时为全实数区间。各自按自己的 miss 更新 α。端点及其后续反馈的总影响不等于作者早期 0.832 的因果分解。',
              '', '| 突变窗口 | A 方法 | 覆盖 ± SE | 均宽 | Winkler | 无穷点/总点 | 含无穷种子/100 | α=0 | α<0.01 | 半宽=历史最大 |',
              '|---|---|---:|---:|---:|---|---|---:|---:|---:|'])
for w in ('early50','transition20','pure30','high'):
    for m in ('A','A_full','A_clip','A_order','A_inf'):
        d=s(m,w);a=r['adaptive_diagnostics'][f'main/shift/{m}/{w}']
        lines.append(f"| {w} | {m} | {stat(d['coverage'])} | {f(d['width']['mean'])} | {f(d['winkler']['mean'])} | {a['infinite_points']}/{a['point_denominator']} | {a['infinite_seeds']}/100 | {a['alpha_zero']} | {a['alpha_lt_001']} | {a['halfwidth_equals_history_max']} |")
lines.extend(['','| 重点端点比较（突变条件） | 窗口 | 覆盖配对差 ± SE | 95% 区间 |','|---|---|---:|---|'])
for name in ('A_inf-A','A_order-A'):
    for w in ('early50','transition20','pure30','high'):
        d=cmp(name,w)['coverage'];lines.append(f"| {name} | {w} | {stat(d)} | {ci(d)} |")
lines.extend(['','无穷半宽按全实数区间记覆盖，均宽/Winkler 为正无穷。只要任一种子窗口含无穷，该比较的普通代价配对差、SE、bootstrap 与偏差比都不可估计；没有删除这些点，没有计算 ∞−∞。两侧有限时才汇总代价差。两条件、全部六窗口的分母及计数均保留于 JSON。',
              '',f"G+L 路径包含 G，所有种子窗口和相同 bootstrap ID 的嵌套单调性检查通过。M-L 原始门槛={r['raw']['M_L']}，仅为杠杆对齐诊断，不是第二次确认。rolling SCP 固定为 100 个绝对残差的第 96 顺序统计量。",
              ''])
fdiag=json.loads((ROOT/'fit_diagnostics.json').read_text());mainfit=[d for d in fdiag if d['setup']=='main']
lines.extend([f"B-df：主设置 df 范围 [{min(d['df'] for d in mainfit):.5f}, {max(d['df'] for d in mainfit):.5f}]，尺度倍数 [{min(d['scale_ratio'] for d in mainfit):.6f}, {max(d['scale_ratio'] for d in mainfit):.6f}]（理论上界 40/37）。G−B-df 原始 M 门槛={r['raw']['M_against_B_df']}。该分母不是无偏尺度估计，不用于替换主 M 或归因作者欠覆盖；全部拟合 df/尺度在 [拟合诊断](evaluation/fit_diagnostics.json)。",
              '', '## 失败、未解原因与范围',''])
failed=[d for d in r['N_F']['items'] if d['state'] in ('不相容','不可估计')]
if failed:
    lines.append('质量门失败项：'+'；'.join(d['condition']+'/'+d['window']+'/'+d['method'] for d in failed)+'。只做静态审计，没有为过门改参数、替换种子、延长窗口或选择 L2。')
    lines.append('独立代码路径复核了观测映射、231 列特征、nλ=16、截距惩罚、训练 RMS、校准独立性、杠杆、测试前残差历史和全部区间；未发现能解释差距的共有实现错误。作者代码/逐重复数据不可得，精确相位、边界索引及 A 约定未核实，无法把统计不相容归因为某个单独原因；仍保留实现未捕获错误的可能。')
else:lines.append('质量门通过仅表示六项诊断未发现预定容差之外的不相容，不证明实现正确或准确标定。')
lines.extend(['原论文已经讨论近期误差和顺序更新。本次不把它描述为首次发现。B/G 除尺度来源外还有杠杆因子、样本内/外残差与噪声水平差异；比较不能隔离在线更新的单独因果效应。A 的独立主约定不匹配原文数值不能证明作者错误。',
              '', '百分位 bootstrap、六项相容区间和固定序列只承诺名义/近似误差控制，不保证有限样本精确覆盖；未追加 BCa、种子或事后阈值。',
              '', '## 复现与交付','',
              '见 [执行说明](EXECUTION.md)、[README](README.md)、[精确依赖](requirements.lock)。本地 `artifacts/v1/` 保存创新、带噪观测、划分、拟合、逐点区间及统一 bootstrap ID；全部原始文件均有哈希且不入库。仓库只含小型结果、种子级汇总、图和元数据；不含论文全文、作者图或大文件。',
              '', '英文 X 帖为 [x_post.md](x_post.md) 的草稿，未发送。仅本地中文提交，不 push。',''])
if audit['status']!='passed':
    lines=[('独立复算仍在运行，尚未出具通过结论；此文为待审核草稿。' if line.startswith('独立 [recompute.py]') else line) for line in lines]
    lines.insert(2,'**草稿：独立复算尚未完成，不是最终交付。**')
Path('REPORT.md').write_text('\n'.join(lines))

posts=[
'Independent NGRC replication, conditional on one frozen trajectory and observation noise. 100 paired noise seeds; primary and all prespecified sensitivity runs completed. This checks a known published gap, not a new discovery.',
('The preregistered baseline/control gate failed, so confirmatory claims are blocked. Raw interval results remain descriptive. A mismatch does not establish an error in the original paper.' if not r['N_F']['passed'] else
 'The preregistered baseline/control gate passed. Main and sensitivity results are reported separately; coverage, closeness to 95%, and interval cost remain distinct questions.'),
'Frozen protocol, code, input/output hashes, interval logs, and independent recomputation are retained. Adaptive CP conventions are uncertain. Prespecified infinite-width variants retain infinite costs; no bad seeds or points were dropped.',
'Limits: inference concerns observation noise on the frozen trajectory, not arbitrary trajectories or real-world guarantees. Exact first-public timestamp remains unverified. Report includes failures, costs, and uncertainty; draft only.'
]
assert len(posts)<=4 and all(len(p)<=270 for p in posts)
Path('x_post.md').write_text('# English X drafts (not posted)\n\n'+'\n\n'.join(f'<!-- Post {i+1}: {len(p)} characters -->\n{p}' for i,p in enumerate(posts))+'\n')
write_json(EVAL/'delivery_metadata.json',{'raw_bytes':sum(d['bytes'] for d in json.loads((ROOT/'artifact_manifest.json').read_text()).values()),
           'raw_file_count':len(json.loads((ROOT/'artifact_manifest.json').read_text())),
           'x_post_characters':[len(p) for p in posts], 'prereg_unchanged':True,
           'first_public_timestamp_verified':False,'pushed':False})
print('Rendered REPORT.md, 3 figures, all diagnostics, JSON, and four length-checked draft posts.')
