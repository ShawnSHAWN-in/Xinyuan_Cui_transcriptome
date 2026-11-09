# eggNOG非模式生物注释背景

从蛋白FASTA的eggNOG-mapper结果构建基因层面的GO/KEGG关联，复用现有PublicDataManager、TERM2GENE/TERM2NAME及富集接口。
本步骤不运行eggNOG、不载入表达矩阵，也不复制Seurat对象。合成测试通过不代表实际物种注释覆盖率已验证。

## 输入与ID

- 上传与定量参考版本一致的蛋白FASTA。网站的真核生物输入应选择蛋白序列。
- 下载文本主表`*.emapper.annotations`或其`.gz`文件；不使用`.hits`、`.seed_orthologs`或Excel替代。
- 支持网站v2/v3共有的列名；读取`#query`或`#query_name`表头，忽略`##`说明及尾部统计，不依赖固定列号。
- 必需列：`query`、`GOs`、`KEGG_Pathway`。可选保留列：`KEGG_ko`、`Preferred_name`、`Description`。
- 每行一个蛋白query；`GOs`和`KEGG_Pathway`为逗号分隔列表，空值、`-`、`NA`表示无注释。重复query及不齐整TSV会停止。
- 解析程序沿用eggNOG结果中已完成的注释转移设置，不额外猜测或更改比对阈值、分类范围及置信度门槛；原输出头记录到metadata。

必须提供蛋白到分析基因的映射表，TSV字段如下：

| 原样列名 | 类型 | 含义 |
| --- | --- | --- |
| query | 字符串 | FASTA序列首个标识符，与eggNOG query完全一致 |
| Geneid | 字符串 | 下游表达矩阵主键；单细胞项目需对应Seurat实际特征名 |

每个query只能对应一个Geneid；多个query可以对应同一基因，注释取并集并去重。建议映射表覆盖上传的完整蛋白集，包含无命中蛋白。
不自动删除`gene-`、蛋白版本号或其他前缀；不使用Preferred_name猜测主键。多个不同基因拥有相同Preferred_name不会被合并。
只有query本来就是分析基因ID时才可明确设置`query_is_gene_id=TRUE`，此时不提供映射表。
这种模式无法从结果表推断未输出蛋白的存在，覆盖率分母只包含结果表中的唯一query；需要完整分母时仍应提供完整映射表。

KEGG名称可提供本地TSV（`TERM`、`NAME`两列字符串），或由KEGGREST在线获取参考通路名称。
`TERM`接受`map00010`、`ko00010`或带`path:`前缀的对应ID，统一为`map00010`；同一通路必须有唯一名称。
程序直接使用eggNOG的KEGG_Pathway关联，不把KO当作基因，也不会使用猫、人等物种的基因背景代替。
本地名称缺失/冲突或网络失败会明确停止；不会编造通路名称。生成的KEGG_TERM2NAME.tsv可作为后续离线名称快照。

## 独立准备参考文件

在仓库根目录逐段运行，先替换真实输入及物种：

```r
source("analysis_scripts/lib/eggnog_annotation.R")
eggnog_file     <- "00.database/eggnog/reference.emapper.annotations"
query_gene_file <- "00.database/eggnog/protein_gene.tsv"
species        <- "填写实际参考物种"
output_dir     <- "00.database/eggnog/prepared"
kegg_names     <- NULL  # NULL首次在线获取；离线时提供TERM/NAME两列TSV
reference <- build_eggnog_annotation_lst(eggnog_file, query_gene_file, species = species,
  ontology = "BP", include_ancestors = TRUE, kegg_term2name_file = kegg_names)
save_eggnog_reference_paths(reference, output_dir)
reference$coverage
rm(reference)
gc()
```

依赖：data.table、tibble、dplyr >= 1.1.1、AnnotationDbi、GO.db；在线名称获取另需KEGGREST。
GO.db用于筛选BP/CC/MF及可选祖先展开，沿用仓库既有算法；未被当前GO.db识别的ID单独输出，需核对注释数据库版本。
基础函数支持BP、CC、MF；当前bulk项目和图像约定只启用BP。

## 既有类接口

在按现有配置载入基座后，原`load_combined_gtf(eggnog=...)`入口现已可用：

```r
manager <- PublicDataManager$new(kegg_organism_code = NULL, cache_dir = "00.database/eggnog/cache")
manager$load_combined_gtf(eggnog = eggnog_file, query_gene_file = query_gene_file, species = species,
  include_ancestors = TRUE, kegg_term2name_file = kegg_names)
manager$get_gene_annotation()
manager$get_go_background("BP")
manager$get_kegg_background()
manager$get_annotation_metadata()
```

项目入口`load_project_annotation(eggnog=..., query_gene_file=..., species=...)`同样支持。
不再要求有效的物种KEGG代码，也不执行物种KEGG联网预检。GTF不能代替蛋白到Geneid映射；此模式不解析或猜测GTF关联。
同次调用不混用OrgDb、gene2go或gene2ensembl。基础类单独使用时需先source本页的eggnog_annotation.R。

为了兼容既有富集和差异结果接口，`ENTREZID`列在eggNOG模式存放**原始Geneid的字符串别名**，不是NCBI GeneID。
metadata明确记录`id_type=Geneid`和`compatibility_id_column=ENTREZID`。不得把这列送入Entrez专用的ID转换或Pathview。
targets的Pathview分支会输出SKIPPED原因；交互配置也应保持`run_pathview=FALSE`。

## 接入已有分析入口

交互脚本：在`analysis_scripts/00_config.R`的cfg中设置以下字段，再运行原01模块：

```r
cfg$annotation_source <- "eggnog"
cfg$eggnog_file <- eggnog_file
cfg$query_gene_file <- query_gene_file
cfg$query_is_gene_id <- FALSE
cfg$kegg_names_file <- kegg_names
cfg$reference_species <- species
cfg$run_pathview <- FALSE
```

保留原实验设计和统计流程；01模块此时只检查eggNOG相关注释输入，不要求NCBI GTF或gene2go文件。
注释输入或配置变化后需重跑01及其依赖模块，不能混用旧背景。

已有targets工程：`config/reference.yml`设置`annotation_source: eggnog`、真实`reference_species`和`query_is_gene_id: false`；
`config/paths.yml`填写`eggnog_file`、`query_gene_file`，可选填写`kegg_names_file`。`run_kegg`仍控制KEGG，默认配置仍为ncbi。
沿用现有annotation_go、annotation_files、kegg_snapshot_files、background_jobs端口，没有另建执行工程。
蛋白映射、名称表及注释文件均由file目标跟踪；不使用文件名猜测版本。

## 输出与字段字典

独立导出位于指定output_dir；现有bulk入口位于`02_annotation/eggNOG/`。
每个TSV都输出表头；空结果仍保留字段。比例为0–1，缺失映射使用NA。

| 文件 | 每行对象 | 全部字段（类型和含义） |
| --- | --- | --- |
| gene_annotation.tsv | 一个映射表中的基因 | Geneid：字符串，分析主键；SYMBOL：字符串，本模式沿用Geneid；ENTREZID：字符串，Geneid兼容别名；GENENAME：字符串，沿用Geneid，不伪造全名 |
| gene2go_direct.tsv | 一条基因–GO直接关联 | ENTREZID：字符串，Geneid别名；GO：字符串，GO编号；SYMBOL：字符串，Geneid |
| query_mapping.tsv | 一个蛋白query | query：字符串，蛋白ID；Geneid：字符串/NA，基因ID，NA表示缺映射；HasAnnotationRow：逻辑，是否出现在结果表；HasGO：逻辑，原GOs非缺失，不代表属于选定本体；HasKEGG：逻辑，原KEGG_Pathway非缺失 |
| coverage.tsv | 一个数据库 | Database：字符串，GO_BP/GO_CC/GO_MF或KEGG；Annotated：整数，有所选数据库关联的唯一基因数；Total：整数，完整映射表中的唯一Geneid数，身份映射时为结果query数；Proportion：数值，Annotated/Total |
| unknown_go.tsv | 一个未知GO编号 | GO：字符串，未在当前GO.db中识别的原始GO编号 |
| gene_map.tsv | 一个基因 | gene：字符串，分析主键；enrichment_id：字符串，同一主键；满足单细胞一对一映射端口 |
| GO_TERM2GENE.tsv | 一条选定本体的基因–GO关联 | TERM：字符串，GO编号；GENE：字符串，Geneid；按include_ancestors决定是否含祖先关联 |
| GO_TERM2NAME.tsv | 一个选定本体的GO条目 | TERM：字符串，GO编号；NAME：字符串，GO.db条目名称 |
| KEGG_TERM2GENE.tsv | 一条基因–通路关联 | TERM：字符串，统一为map编号；GENE：字符串，Geneid；多个蛋白和map/ko重复已去除 |
| KEGG_TERM2NAME.tsv | 一个实际关联的通路 | TERM：字符串，map编号；NAME：字符串，来源于本地快照或KEGGREST的名称 |
| gene_sets.rds | 数据库具名list | GO_BP（或指定本体）、KEGG；无关联数据库省略；每个成员含term2gene(term_id,gene_id)、term2name(term_id,term_name)、species、id_type、source；列均字符串，元数据为单个字符串 |
| annotation_metadata.rds | 注释来源list | source：注释工具；species：参考物种；id_type：实际基因ID类型；compatibility_id_column：兼容列名；ontology：本体；include_ancestors/run_kegg：逻辑参数；header：原始说明行；files：输入绝对路径命名的MD5；go_db_version：GO.db版本；kegg_names：名称来源；created：生成时间 |

映射表中无注释命中的蛋白仍保留在query_mapping；无注释的基因不填假GO或假通路。
注释背景是功能关联库，不等于每次差异比较的统计universe。ORA仍应将当前比较的全部受检/可检测基因与对应背景取交集，不能只用显著DEG当分母。
上传全部参考蛋白，避免仅给差异基因注释；多个蛋白合并后按基因计数。

## 验证与来源

```bash
Rscript analysis_scripts/tests/eggnog_annotation.R
```

合成测试覆盖v2/v3列顺序、gzip、注释尾部、显式映射、蛋白合并、BP筛选/祖先展开、KEGG去重、错误输入、输出端口及既有类接口。
测试不运行实际样本差异分析，不证明真实小熊猫的注释覆盖情况；已有完整bulk回归套件需在项目环境和既有fixture齐备后执行。

本次验证：R 4.3.3、dplyr 1.1.4、GO.db 3.18.0下7组合成测试通过；30个相关R文件解析通过。
额外读取官方v2 `test_mmseqs.emapper.annotations`（2个query，blob fb0746616ddbfa4630f86797e564d042ea65acb9）及
v3 `test_diamond.emapper.annotations`（20个query，blob 276df880620d21c00473a85896e4dd46aa108311）成功。
官方fixture仅用于格式读取检查，未将其中的生物结果作为本项目背景；无表头文件明确拒绝，不猜列号。
类接口测试单独载入实际类定义和兼容方法，未执行整个bulk流程；未执行targets完整调度、服务器实际注释或真实富集分析。

- [网站提交页面](https://eggnog-mapper.cgmlab.org/submit/)
- [官方v3输出说明](https://github.com/eggnogdb/eggnog-mapper/blob/main/USAGE.md#output-files)
- [官方v2输出说明](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.13)
