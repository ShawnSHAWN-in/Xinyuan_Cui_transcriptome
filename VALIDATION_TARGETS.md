# targets 工程验证记录

验证使用仓库脚本生成的合成数据，400个基因、36个样本、九个默认比较；不包含用户真实表达矩阵。工程化代码复用原R6 edgeR计算函数。

## 2026-10-20 绘图恢复验收

- 元数据、原有绘图合同和新增旧版绘图回归测试通过。
- 火山图、MA、rank逐项对照原R6绘图函数：底层数据、标签层、绘图参数和颜色一致，标签数沿用8。
- 实际导出九个比较的27组差异单图和3组差异组合图，各有PDF/PNG。
- QC对照旧模块验证filterByExpr分组过滤、未缩放PCA、相关性及归一化因子；元数据行顺序改变后仍正确按样本对齐。
- 默认QC实际导出33组PDF/PNG（含5组兼容副本，3D包可用），并检查文件存在、非空、PNG可解码且不超过像素预算。
- 差异图和QC共126个PDF/PNG全部非空，63个PNG均成功解码；核验PNG结束块及PDF结束标记。
- 导出先写临时文件，通过完整性检查后再替换正式文件；测试确认绘图失败或临时文件损坏时保留已有完整图像。
- 查看火山图、MA、rank、四模板PCA组合、样本概览及3D得分图，核对标签、颜色、图例和布局。
- edgeR拟合、注释关联、DEG筛选、GSEA排序四个统计函数的函数体和参数与上一工程提交相同。

本次局部验证环境为R 4.3.3、edgeR 4.0.16、ggplot2 3.4.4及ggpubr 0.6.0；保留原renv.lock。本次针对绘图恢复进行验证，没有重新运行联网KEGG/Pathview或用户服务器真实数据。

```bash
Rscript tests/generate_baselines.R --inputs-only
Rscript tests/smoke_restored_plots.R
```

本次图像和配置说明见 [docs/QC_OUTPUTS.md](docs/QC_OUTPUTS.md)。以下321任务等记录属于2026-10-19工程首版，保留用于数值基线追踪。

## 2026-10-19 工程首版验收


- 单元验证全部通过，包括九比较数值基线、元数据分组与绘图合同。
- 原九比较的受检基因顺序、logFC、P值、FDR、logCPM、baseMean与工程化前保存的基线一致。
- 默认完整任务图完成321个执行任务，覆盖9个edgeR、36个ORA和18个GSEA，以及QC、表达热图、集合图、富集单图和组合图。
- 与旧流程富集结果逐ID对齐：ORA P值/FDR最大绝对差约2.5×10^-32；GSEA NES最大绝对差约4.9×10^-15，P值/FDR最大绝对差约5.0×10^-16。
- 原样再次执行全部跳过。动态聚合节点和全局函数的元数据更新时间不能用来判定任务是否重算，测试使用任务进度及具体统计分支时间。
- 通用分类变量/分类时间组合、factor水平顺序、分层与显式组对、缺失对照报错、Up/Down固定位置、数据库不可用标注和PNG像素预算已进行单元验证。
- 查看合成PCA和富集比较图，确认白底黑框、方向形状、图例、长标签及时间分隔。预览仅展示风格，不代表实验结果。

- 修改主题后统计分支时间不变；恢复主题再单独删除一个PNG，只完成1个对应绘图分支，edgeR/ORA/GSEA缓存均未改动。
- 将DEG阈值改为|log2FC| >= 3后重跑ORA，上游edgeR拟合与完整排序GSEA分支均复用。
- 启用time_within_condition后得到17个比较，只新增8个时间比较拟合；原9个拟合的任务名与完成时间保持不变。
- 增量冒烟测试和新增分支测试均以PASS退出。

合成输入与回归结果表不提交到仓库。生成脚本直接调用旧模块01/03/06，跳过绘图，不调用新工程统计函数；生成结果已与先前独立保存的基线核对，edgeR逐值相同，GSEA NES最大绝对差约1.0×10^-15。

## 环境

R 4.3.3 / Bioconductor 3.18；targets 1.12.0、edgeR 4.0.16、clusterProfiler 4.10.1、fgsea 1.28.0、ggplot2 3.5.2。完整依赖快照见renv.lock。

本地测试执行 `tar_make(callr_function = NULL)`。测试运行环境的PID命名空间与/proc不一致，进程监控库在本地验证运行时按/proc/self确定当前监控PID；此处理没有写入工程代码或上传的包源码。服务器使用标准R环境，不需要该兼容处理。

本次KEGG使用带表头的本地合成TERM2GENE/TERM2NAME快照。实时KEGG下载、Pathview联网绘制及用户服务器上的真实数据执行未在此验证。

## 重复验证

```bash
Rscript tests/generate_baselines.R
Rscript tests/run_tests.R
Rscript tests/smoke_targets.R
TARGETS_TEST_ROOT=/path/to/smoke/TEST_ROOT Rscript tests/check_branching.R
```

完整冒烟测试默认包含网络图和UpSet；单条目UpSet可能产生geom_line只有一个观测的提示，图形仍可完成。

