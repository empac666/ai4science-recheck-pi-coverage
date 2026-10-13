# NEXT：NGRC 预测区间独立复核

本仓库已完成 NEXTrun 的预定主实验与全部 L2；正式结果与研究结论见 [中文报告](REPORT.md) 和 [完整结果 JSON](evaluation/result.json)。结果是否支持命题由锁定质量门及 M→E→J 顺序决定，完成运行不等于假设成立。

- 锁定协议：[PREREG.md](PREREG.md)，原件不改；文首“尚未公开”是锁定时的历史状态。
- 结果前执行冻结：[EXECUTION.md](EXECUTION.md)、[源文件哈希](EXECUTION_SHA256SUMS)、[配置](config.json)、[软件锁](software_lock.json)。
- 验收与复算：[结果前验收](validation_before_run.txt)、[独立复算脚本](scripts/recompute.py)、[复算记录](evaluation/independent_recompute.json)。
- 审计：[偏离记录](DEVIATIONS.md)、[公开读取证据](evidence/public_lock.json)、[300 次尝试](evaluation/attempts.json)、[原始文件哈希](evaluation/artifact_manifest.json)。
- 交流：[中文 X 单帖草稿](x_post_cn.md)、[英文 X 串帖草稿](x_post.md)，均未发送；全部为本地中文提交，不 push。

本研究只在冻结轨道 x* 的观测噪声总体内推断，是对已知原表差距的独立实现检验。作者代码与逐重复数据不可得；A 使用本协议主约定，不能因一项不匹配就认定原论文错误。首次公开的精确时间未能从 GitHub events 接口核实；本次保存了结果生成前的匿名公开可读证据，不用本地 Git 时间或本次读取时间冒充首次公开时间。

## 复现执行

精确执行提交为 `c061a7eafaa1073e0e07bd1a5da3aae3f2fbae62`，这是首次正式数据生成前的实现冻结。正式执行器要求工作树干净、依赖版本和冻结文件哈希一致，且拒绝覆盖原目录。可在另一个目录检出该提交，安装 [requirements.lock](requirements.lock) 中的依赖；Python/BLAS 信息见软件锁。当前工作目录的 `.venv` 复用了固定系统包，限制已登记。

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 7200 .venv/bin/python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 7200 .venv/bin/python -m nextrun.run --output artifacts/v1
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 7200 .venv/bin/python scripts/recompute.py artifacts/v1 --output artifacts/recompute_v1.json
```

上述命令在新检出的执行冻结目录运行。报告渲染器随后加入交付提交，只做展示，不改变封存统计；在含该脚本的交付目录可运行：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 7200 .venv/bin/python scripts/render_report.py
```

`artifacts/v1/` 原始数据只留本地，所有文件都有 SHA-256；它包含标准正态创新、带噪观测、划分、拟合、逐点区间、预测前 α/历史极值以及统一 bootstrap ID。原始大文件不入库。仓库中的 `evaluation/seed_metrics.npz` 为小型种子级汇总，原始数组可以由冻结代码生成；报告的独立复算需要原始本地产物或按冻结代码重建的完整目录。

`SHA256SUMS` 校验最终交付文件；`EXECUTION_SHA256SUMS` 校验结果前执行冻结。后者如涉及追加式审计文档，应在记录的执行提交上校验，执行核心代码从冻结后保持不变。

## 历史锁定流程

下述流程是 NEXTlock 阶段保留的历史说明；“本轮”指锁定阶段，执行状态以上述 NEXTrun 结果为准。

## 锁定流程

1. 本轮只做文档静态检查、更新 `SHA256SUMS` 并创建中文本地提交。读取完整提交 SHA 和 `PREREG.md` 的 SHA-256；状态保持“本地定稿、待公开锁定”。本地 commit 时间、文件修改时间和预先书写的日期均不构成公开预注册时间戳。
2. 后续获授权发布时，**以本轮定稿提交作为首次 push 到公开 GitHub 的分支头；该首次公开 push 的提交是预注册锁定提交，首次公开可见时间是时间戳**。发布前确认公开仓库可访问且该提交仍无正式实验结果，核对协议哈希。若发布前确需改协议，先形成有说明的新定稿提交，再将该新提交作为首次公开 push 的锁定提交；不能用后来的提交冒充此前公开版本。本任务不执行 push。
3. 首次公开后，用单独的后续登记提交保存完整锁定 SHA、对应 `PREREG.md` SHA-256、公开提交固定链接、GitHub 可核验的首次公开时间及其证据链接。只记录必要的公开版本信息，不收录个人或主机标识、凭据或本地路径。不用 Git author/committer 日期代替首次公开时间，不以私有仓库提交时间回溯公开登记；若无法核实公开时间，明确标为待核实。登记提交只引用锁定提交，不改写该快照，也不把自身当成第一次锁定。
4. 公开锁定之后、正式数据生成之前，完成独立实现、协议规定的非实验验收、配置及精确软件依赖冻结，并另行提交执行版本，引用锁定 SHA/协议哈希。100 个正式种子、模拟数据、拟合与结果计算须等公开锁定和执行冻结均完成后才开始；本轮没有进行这些工作。
5. 公开后任何协议改动均另建修订提交，写明原因、改动范围及是否已接触结果，保留原锁定提交/哈希和偏离记录；结果后修复按修订分析报告，不改写历史、替换时间戳或伪装成原预注册分析。

协议内不嵌入其自身 SHA-256 或包含自身内容的提交 SHA，避免循环引用；最终值从 Git 和文件读取，随后在独立登记记录中引用。未来 0.5–2 小时计算预算只是静态估计，未做性能试跑。

论文及作者材料保留原权利；仓库只分发原创复核文档、来源元数据与链接，不打包 PDF、全文 HTML、TeX、原图或作者代码，不替第三方材料授予开放许可。

## English

The preregistered primary experiment and all scheduled sensitivity analyses have been executed. See REPORT.md and evaluation/result.json for the actual quality gate and confirmatory decisions. Completion does not imply confirmation. Inference is conditional on observation noise around the frozen trajectory; this is an independent implementation check of a known published gap.

The original PREREG.md is unchanged. The implementation, tests, configuration and precise dependencies were committed before formal seed generation. Raw input/output files remain outside Git with checksums; a standalone verifier re-integrates the trajectories and reconstructs interval metrics and statistical decisions. All commits are local; no push or social post was made. Exact first-public time is unverified; only pre-execution anonymous public readability is established.

Third-party paper texts and original figures are not distributed. The repository contains our code, analysis, redrawn figures, metadata and source links.
