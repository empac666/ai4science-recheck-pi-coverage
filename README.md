# NEXT：NGRC 预测区间独立复核 / Independent NGRC interval replication

**状态：结果前协议；结果未出；实验尚未运行。当前仅为本地仓库雏形，未推送，未取得公开预注册时间戳。**

本项目计划对 [arXiv:2609.32169v1](https://arxiv.org/abs/2609.32169v1)，*Uncertainty Quantification of Next Generation Reservoir Computing with Applications to Memory-Driven Dynamical Systems*，进行**预注册独立复核**。范围仅限 §5.4 的合成观测噪声突变实验：比较冻结 Bayesian 区间、冻结 split conformal、更新 Gaussian 与 adaptive conformal 的覆盖、宽度和 Winkler 评分。

这里冻结的是结果前研究协议，尚无实验实现、数据或复核结果。文档中的 0.726、0.929、0.938 等数字来自原论文表 1，均不是本项目的测量值。本项目独立于论文作者；作者代码与逐种子数据尚未核实获得，不承诺逐位复现。

计划以 100 个配对观测噪声重复评价同一条固定潜在轨道。全部统计推断都条件于该轨道及预定设置，不能推广为任意轨道或工业场景的覆盖保证。主收益判据要求 G−B 的双侧 95% 配对区间下界大于 0.15；确认性检验按 M→E→J 固定顺序进行，负对照失败会阻断确认性结论。

| 文件 | 内容 |
|---|---|
| [PREREG.md](PREREG.md) | 冻结的观测映射、索引、方法、统计判据、敏感性与失败规则 |
| [feasibility.md](feasibility.md) | 来源核查、许可范围、依赖与静态计算预算 |
| [response.md](response.md) | N-C01 至 N-C11 的逐条审查回复 |
| [evidence/README.md](evidence/README.md) | 随包来源元数据及未分发全文的边界 |
| [SHA256SUMS](SHA256SUMS) | 本仓库文档与元数据的 SHA-256 清单 |

下一阶段须先编写独立实现，完成协议规定的验收并冻结配置、精确依赖和执行版本，再开始正式实验。当前 0.5–2 小时计算预算只是静态估计，未做性能试跑。任何协议变更应留下新提交和说明，不覆盖本快照。

论文及作者材料的权利归原权利人；本仓库只提供原文链接、来源元数据与原创复核文档，不打包 PDF、全文 HTML、TeX、原图或作者代码，也不替第三方材料授予开放许可。

## English

**Status: pre-results protocol. No experiments have been run and no replication results exist. This is a local repository scaffold; it has not been pushed or publicly preregistered.**

This project is a **planned preregistered, independent replication** of the synthetic observation-noise shift experiment in §5.4 of [arXiv:2609.32169v1](https://arxiv.org/abs/2609.32169v1), *Uncertainty Quantification of Next Generation Reservoir Computing with Applications to Memory-Driven Dynamical Systems*. It compares frozen Bayesian, frozen split-conformal, updated Gaussian, and adaptive conformal intervals using coverage, width, and the Winkler score.

The repository freezes a protocol before new results. It contains no experimental implementation, generated data, or replication findings. Numbers such as 0.726, 0.929, and 0.938 are the original paper's reported values, not our measurements. This work is independent of the paper's authors. Author code and per-replication data have not been verified as available, so bitwise reproduction is not promised.

All inference is conditional on one fixed latent trajectory and the specified setup, with uncertainty across 100 paired observation-noise replications. It provides no guarantee over arbitrary trajectories or industrial settings. The primary improvement claim requires the lower endpoint of a two-sided 95% paired interval for G−B to exceed 0.15. Confirmatory testing follows the fixed sequence M→E→J; failed negative controls block confirmatory claims. Secondary diagnostics cannot replace a failed primary comparison.

Read [PREREG.md](PREREG.md) for the authoritative protocol, [feasibility.md](feasibility.md) for scope and budget, and [response.md](response.md) for the review responses. These technical documents are in Chinese. [Evidence metadata](evidence/README.md) and [SHA256SUMS](SHA256SUMS) support provenance and integrity checks. A local Git commit and file hash do not constitute a third-party public preregistration timestamp.

Before any formal experiment, the implementation, specified validation checks, configuration, exact dependencies, and execution commit must be frozen. Runtime estimates have not been benchmarked. Protocol changes must be documented in new commits. The paper and third-party materials retain their original rights; this repository does not redistribute their full text, figures, or code.

## 文件完整性 / Integrity

以下命令只校验文件或读取 Git 信息，不运行实验。These commands only inspect files and Git metadata; they do not run experiments.

```sh
sha256sum -c SHA256SUMS
sha256sum PREREG.md
git rev-parse HEAD
```
