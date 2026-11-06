# 合成测试夹具

本目录全部为软件回归测试数据，不包含用户实验计数或实验结果。所有输入表和结果表均在本地生成，不提交到GitHub。首次测试前运行 `Rscript tests/generate_baselines.R`，直接调用旧模块01/03/06生成合成输入及本地基线；此脚本跳过绘图，不调用新工程统计函数。

- `counts.tsv`：400个虚构Gene0001–Gene0400 × 36个模拟样本，Length均为1000。生成代码在 `analysis_scripts/tests/smoke_core.R`：`set.seed(21)`，先用runif构造均值并施加已知效应，再用rnbinom生成整数计数。
- `samples.tsv`：与模拟矩阵匹配的实验编码模板。
- `reference.gtf`：chr1上的规则位置，虚构Geneid及9000001起的模拟GeneID；不是NCBI真实参考序列的节选。
- `gene2go.tsv`、`KEGG_TERM2GENE.tsv`、`KEGG_TERM2NAME.tsv`：将模拟GeneID分配给少量有效GO/KEGG条目，仅用于覆盖软件分支，不能用于生物学解释。
- `legacy_edgeR_expected.tsv`：原模块03对上述合成计数计算的9个比较、3600行基线，用于检测迁移是否改变统计结果。
- `legacy_enrichment_expected.tsv`：旧富集模块对上述模拟映射计算的ORA/GSEA基线。

已逐元素重新生成并核对全部14,400个计数；edgeR基线重新从旧合成冒烟结果导出，与此处文件字节一致。单元测试继续把新工程结果与这些独立保存的基线比较。
