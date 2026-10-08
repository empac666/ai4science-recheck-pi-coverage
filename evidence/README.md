# 来源元数据 / Source metadata

这些记录来自 2026-10-08 的前轮匿名来源核查。它们记录下载与文件清单，不含新实验结果；NEXTfix 没有重新下载论文、生成模拟数据或拟合模型。

- [probe_manifest.json](probe_manifest.json)：固定论文版本 URL、下载时间、HTTP 状态、字节数、SHA-256，以及当时读取的包元数据。许可页首次 404 如实保留。
- [probe_license_manifest.json](probe_license_manifest.json)：通过官方页面纠正许可地址后的下载记录。
- [source_inventory.json](source_inventory.json)：源包条目清单；该包包含论文素材，不能据此称已获得可执行实验代码。

上述三个 JSON 按原样复制。它们描述的 PDF、HTML、TeX 源包、检索响应或其他全文副本**不在本仓库内**，只保存在上游工作目录供离线核查。字段中的文件名是历史记录，不是声明这些文件随仓库分发。公开读者可由清单中的固定版本 URL 访问来源；下载字节变化时须报告新哈希，不能静默替换协议记录。

观测映射 `y_i=x(t_i)+ε_i` 核对自论文正文 §5.4；三状态 Euler 递推核对自附录 C.5.2。论文链接及关键文件哈希另见 [PREREG.md](../PREREG.md)。[仓库校验清单](../SHA256SUMS)覆盖实际随包文件，不声称包含全部本地证据。

These unchanged metadata files document the earlier anonymous source audit. The full paper, HTML, TeX archive, figures, search responses, and other third-party content are not distributed here. Filenames inside the metadata identify historical downloads, not bundled assets. This repository contains no new experimental results. Original source URLs and hashes are retained for provenance; local Git history is not a public preregistration timestamp.
