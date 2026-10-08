# 执行冻结

协议父提交：`ed8d08bdd29e6084ccb42345896352ea1e62ad8e`；协议 SHA-256：`9d486db80eb5faca0d5aeb21c60ae80088991514f2360fd99208b22e71c77991`。`PREREG.md` 保持原字节。

本文件及 `EXECUTION_SHA256SUMS` 所列代码、测试、配置、依赖锁在正式数据生成前提交。实际执行提交 SHA 由运行器读取 Git HEAD 写入 `run_record.json`，避免自包含哈希循环。主执行器拒绝脏工作树、冻结文件哈希不符或覆盖已有输出目录。所有正式种子首次运行后封存原始输出；独立复算器不导入 `nextrun`，仅从原始观测、拟合与逐点区间日志重算，输出另存。

- RNG：NumPy `Generator(PCG64)`，创新与 bootstrap 流严格依协议；无非预定正式种子。
- 精度：float64；所有固定配置在 `config.json`，实现中的协议常数由验收检查覆盖。
- 单进程、BLAS/OMP/MKL/NumExpr 均 1 线程，无 GPU；正式执行内存地址空间上限 1 GiB，产物预算 1 GB。
- 源码使用 NumPy Cholesky 和线性方程求解；有效自由度求解 `L Z=I`，不构造显式逆或 1600×1600 帽子矩阵。
- 程序先完成主设置判定，再计算敏感旗标；不因质量门失败删减预定 L2。
- 非正式验收仅用人工数组与种子 17 的人工回归设计，未使用 100 个正式种子、正式轨道或试跑结果调参。
- 测试记录见 `validation_before_run.txt`；软件锁见 `software_lock.json`、`requirements.lock`；公开证据见 `evidence/public_lock.json`；限制见 `DEVIATIONS.md`。

执行与独立复算命令均设外部两小时硬限，原始数据只保存于忽略目录，不入库：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 7200 .venv/bin/python -m nextrun.run --output artifacts/v1
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 7200 .venv/bin/python scripts/recompute.py artifacts/v1 --output artifacts/recompute_v1.json
```

中断后目录可能保留 `running` 状态，必须按未完成处理；不得删目录重跑。失败后的修订使用新目录并登记。终止信号由外部硬限处理，成功子集不得替代完整 100 次统计。
