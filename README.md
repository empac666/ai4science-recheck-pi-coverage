# NEXT：NGRC 预测区间独立复核 / Independent NGRC interval replication

**状态：NEXTlock 结果前协议已本地定稿；实验未运行，结果未出，未 push，尚无公开预注册时间戳。**

本项目计划对 [arXiv:2609.32169v1](https://arxiv.org/abs/2609.32169v1)，*Uncertainty Quantification of Next Generation Reservoir Computing with Applications to Memory-Driven Dynamical Systems*，进行预注册独立复核。范围仅限 §5.4 的合成观测噪声突变实验，比较冻结 Bayesian（B）、冻结 split conformal（F）、更新 Gaussian（G）与 adaptive conformal（A）的覆盖、宽度和 Winkler 评分。

文档中的 0.726、0.929、0.938 等是原论文表 1 的数字，不是本项目测量值。作者代码与逐种子数据尚未核实获得，不承诺逐位复现。A 采用本稿预先固定的主约定。本项目是对已知原表差距的独立实现检验，不是盲分析。

全部推断条件于冻结轨道 x* 的观测噪声总体，来自 100 个配对噪声重复；确认性窗口只有测试索引 `{160,…,799}`。M 要求 G−B 的双侧 95% 配对区间下界严格大于 0.15；固定序列为 M→E→J。基准/负对照质量门使用六项 Bonferroni 调整后的统计相容区间，包含 B 的高噪声基准。J 检验两项配对线性差 `G−1.10A` 的 90% 区间上界均小于 0；均值比只作描述。敏感性按协议限定结论措辞，不将变体升级为确认性结果。

| 文件 | 内容 |
|---|---|
| [PREREG.md](PREREG.md) | NEXTlock 最终协议；算法、质量门、统计判据与报告规则的权威版本 |
| [response2.md](response2.md) | 两份二轮复核的 14 项逐条处理；与任务要求的上游同名交付件按字节一致 |
| [response.md](response.md) | 第一轮历史回复；其中被二轮修改的规则以 PREREG 为准 |
| [feasibility.md](feasibility.md) | 来源、许可范围、依赖与静态预算 |
| [evidence/README.md](evidence/README.md) | 来源元数据及全文不随包分发的边界 |
| [SHA256SUMS](SHA256SUMS) | 实际随包文档与元数据的 SHA-256 清单（不自包含） |

## 锁定流程

1. 本轮只做文档静态检查、更新 `SHA256SUMS` 并创建中文本地提交。读取完整提交 SHA 和 `PREREG.md` 的 SHA-256；状态保持“本地定稿、待公开锁定”。本地 commit 时间、文件修改时间和预先书写的日期均不构成公开预注册时间戳。
2. 后续获授权发布时，**以本轮定稿提交作为首次 push 到公开 GitHub 的分支头；该首次公开 push 的提交是预注册锁定提交，首次公开可见时间是时间戳**。发布前确认公开仓库可访问且该提交仍无正式实验结果，核对协议哈希。若发布前确需改协议，先形成有说明的新定稿提交，再将该新提交作为首次公开 push 的锁定提交；不能用后来的提交冒充此前公开版本。本任务不执行 push。
3. 首次公开后，用单独的后续登记提交保存完整锁定 SHA、对应 `PREREG.md` SHA-256、公开提交固定链接、GitHub 可核验的首次公开时间及其证据链接。只记录必要的公开版本信息，不收录个人或主机标识、凭据或本地路径。不用 Git author/committer 日期代替首次公开时间，不以私有仓库提交时间回溯公开登记；若无法核实公开时间，明确标为待核实。登记提交只引用锁定提交，不改写该快照，也不把自身当成第一次锁定。
4. 公开锁定之后、正式数据生成之前，完成独立实现、协议规定的非实验验收、配置及精确软件依赖冻结，并另行提交执行版本，引用锁定 SHA/协议哈希。100 个正式种子、模拟数据、拟合与结果计算须等公开锁定和执行冻结均完成后才开始；本轮没有进行这些工作。
5. 公开后任何协议改动均另建修订提交，写明原因、改动范围及是否已接触结果，保留原锁定提交/哈希和偏离记录；结果后修复按修订分析报告，不改写历史、替换时间戳或伪装成原预注册分析。

协议内不嵌入其自身 SHA-256 或包含自身内容的提交 SHA，避免循环引用；最终值从 Git 和文件读取，随后在独立登记记录中引用。未来 0.5–2 小时计算预算只是静态估计，未做性能试跑。

论文及作者材料保留原权利；仓库只分发原创复核文档、来源元数据与链接，不打包 PDF、全文 HTML、TeX、原图或作者代码，不替第三方材料授予开放许可。

## English

**Status: NEXTlock protocol finalized locally, before results. No experiments have been run, no push has been made, and no public preregistration timestamp exists.**

This is a planned independent replication of the synthetic observation-noise shift experiment in §5.4 of [arXiv:2609.32169v1](https://arxiv.org/abs/2609.32169v1). Published table values are references, not new measurements. Inference concerns observation noise conditional on the frozen trajectory x*, test indices `{160,…,799}`, and the protocol's primary A convention. This is an independent implementation check of an already known published gap.

M requires the lower endpoint of a two-sided 95% paired interval for G−B to exceed 0.15. The confirmatory sequence is M→E→J, preceded by six Bonferroni-adjusted compatibility checks including the B baseline. J requires both paired linear contrasts G−1.10A, for width and Winkler score, to have negative upper endpoints of their 90% intervals. Ratios remain descriptive. Sensitivity findings qualify the reporting language without creating new confirmatory analyses.

The [Chinese locking workflow](#锁定流程) is authoritative: the finalized protocol commit must be the branch tip of the **first push to public GitHub**. That commit anchors preregistration; its first verifiable public availability supplies the timestamp. A local commit date is not a public timestamp. A subsequent registration commit records the lock SHA, PREREG SHA-256, public commit link, public time and evidence without rewriting the locked snapshot. Public locking and a separate implementation/configuration/dependency freeze must precede formal data generation. Amendments remain separate, disclose whether results were seen, and never replace the original timestamp. This task performs no push or experiment.

[PREREG.md](PREREG.md) is the final protocol; [response2.md](response2.md) resolves both second-round reviews. [response.md](response.md) is historical. Full third-party texts, figures and code are not distributed. Runtime estimates have not been benchmarked.

## 文件完整性 / Integrity

以下命令只读取文件和 Git 元数据，不运行实验。These commands only inspect files and Git metadata.

```sh
sha256sum -c SHA256SUMS
sha256sum PREREG.md
git rev-parse HEAD
```
