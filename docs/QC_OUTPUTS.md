# QC和差异图恢复说明

本次按工程化前的 `analysis_scripts/02_sample_qc.R`、`03_differential_expression.R` 和原R6绘图函数恢复。targets构建、灵活分组、富集Up/Down同图及流式导出继续使用。

## 差异图

直接调用原 `run_deg_visualization()`，三种图均沿用旧版8个标签设置。Up为 `#C94C4C`，Down为 `#377EB8`，非显著为 `#C9C9C9`；富集图的配色配置不会改变这三种图。

| 图像 | 标签与颜色规则 | 坐标 |
|---|---|---|
| 火山图 | 显著基因按padj选Top 8，SYMBOL文本标注；按FDR及logFC阈值区分Up/Down/非显著 | logFC与−log10(padj) |
| MA图 | 原ggpubr::ggmaplot；按padj选Top 8；保留黑色粗体矩形标签 | Log2 mean expression与logFC |
| rank图 | 按绝对logFC选Top 8；标签及点按logFC正负着色，与旧版一致 | logFC降序排名与logFC |

单比较三类图和每方案三类组合图均输出PDF/PNG。每个比较方案内使用共同的logFC范围；火山图同时使用共同显著性范围。绘图阈值仍读取analysis.yml，不改变差异分析计算。

## QC计算

- CPM/logCPM沿用edgeR TMM归一化。
- PCA和样本相关性恢复按组合分组执行的 `edgeR::filterByExpr()` 过滤。
- PCA去除零方差特征；默认 `center=TRUE, scale.=FALSE`。
- 四种PCA模板直接复用原基座，使用同一批过滤基因和同一组坐标。
- 默认组合组采用旧版Condition颜色随Time_f加深的梯度。
- 全局分面PCA共用全局坐标；独立分层PCA在每层内重新计算，二者含义在标题中区分。

分类变量和顺序继续读取variables.yml。QC变量在figures.yml单独配置：

```yaml
qc:
  group_by: Condition_Time
  color_by: Condition
  shape_by: Time_f
  strata: [Time_f]
  pca_scale: false
```

可换成其它已声明分类变量或组合变量；strata可以包含多列。数值时间与分类时间分别使用Time_h/Time_f；QC形状和分层使用分类列。factor水平和display.order_by决定展示顺序。

## 单图清单

默认文件位于 `04_targets_analysis/03_sample_qc/`。下表主名均各有PDF和PNG。

| 类别 | 主名 | 默认数量 |
|---|---|---:|
| 文库、检出与归一化 | QC_CountSum、QC_DetectedGenes、QC_ZeroProportion、QC_NormalizationFactor | 4 |
| 表达分布 | expression_boxplot_counts、expression_boxplot_log_cpm、expression_boxplot_tpm、expression_boxplot_fpkm | 4 |
| 四种PCA模板 | PCA_points_ellipse、PCA_sample_labels、PCA_group_centroid、PCA_minimal_points | 4 |
| 全局PCA与分面 | PCA_treatment_time、PCA_global_coordinates_by_time | 2 |
| PCA方差解释 | PCA_explained_variance | 1 |
| 独立分层PCA | PCA_independent_加分层变量与水平 | 默认3 |
| 样本相关性 | sample_correlation_fixed_order、sample_correlation_clustered | 2 |
| 样本聚类 | sample_dendrogram | 1 |
| 3D PCA得分 | PCA_3D_scores | 1，可选依赖 |
| 注释覆盖 | annotation_coverage | 1 |

独立分层PCA同时保存每层坐标及方差解释表。四模板PCA、局部PCA和3D的状态写入 `QC_status.tsv`。原轮子对少于4个样本的组不画椭圆，默认每组合组3个重复时会保留提示；不强行添加不可靠的椭圆。

3D使用真实得分坐标，安装scatterplot3d后自动输出；缺包或主成分不足时记录SKIPPED。KEGG不可用时覆盖率保留NA及数据库状态，避免显示为0%。

## 组合图清单

均位于 `04_targets_analysis/09_combined_figures/`，各有PDF和PNG。

| 主名 | 面板 |
|---|---|
| 01_sample_overview | 四项文库与检出指标 |
| 02_expression_distributions | counts、logCPM、TPM、FPKM四种箱线图 |
| 03_PCA_four_styles | 四种原始PCA模板 |
| 04_independent_PCA_by_time_页码 | 每层独立PCA，较多分层时分页 |
| 05_sample_relationships | 全局PCA、固定顺序相关性热图、聚类树 |

同时重写上一工程版的五组兼容文件：01_PCA、02_library_size、03_expression_boxplot、04_correlation、05_clustering。它们分别对应上述恢复后的图，不是额外的独立分析。

默认36样本、3个时间水平且3D可用时，共有33组PDF/PNG（66个文件，含5组兼容副本）。实际清单以 `QC_figure_manifest.tsv` 为准。矩阵、样本汇总、过滤标志、PCA坐标及方差、相关性和注释覆盖表另存TSV。

## 更新与运行

在工程分支更新代码后运行：

```r
targets::tar_make(names = c(qc_files, deg_plot_files, deg_combined_files))
```

已有缓存有效时复用edgeR和富集统计；本次更换QC过滤会更新QC计算及引用QC数据的表达热图。完整更新可运行 `targets::tar_make()`。

导出继续逐面板绘制并限制PNG像素。每个文件先写入同目录临时路径，关闭设备并检查PNG结束块或PDF结束标记后才替换正式文件；失败时保留已有完整文件。本次修改了共享导出函数，因此完整构建还可能重新导出其它模块的图像。

只验证恢复的绘图模块：

```bash
Rscript tests/generate_baselines.R --inputs-only
Rscript tests/smoke_restored_plots.R
```

测试全部使用本地生成的合成数据，不上传输入矩阵、回归结果表或真实实验结果。
