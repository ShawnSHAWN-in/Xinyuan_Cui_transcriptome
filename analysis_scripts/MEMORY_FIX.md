# 09_combine_export.R 内存修复

本次仅修改模块09；已有01–08模块的分析结果可以继续使用。无需重新运行差异分析、ORA、GSEA或模块07。

## 原因与证据

1. 旧模块09首先读取整个 `07_enrichment_plots.rds`，但实际只需要少量GSEA曲线和Pathview文件路径。该RDS保存所有数据库、比较和图形类型的ggplot、gglist及其关联环境。与当前会话已存在的对象叠加，重复反序列化会增加内存。
2. 本地旧合成数据的绘图RDS约11 MB。预载ggplot2和patchwork后，读取它使RSS从103840 KiB增加到739956 KiB，R堆从41.4 MiB增加到688.8 MiB。`object.size()`仅报告20.44 MiB，不能用它估计反序列化后的完整进程占用。包加载也可能贡献部分增量，以上未进一步拆分。
3. 旧Pathview代码通过 `lapply()` 同时解码整组PNG，默认得到浮点数组，组合grob进一步持有图像。在3000×2000像素RGBA测试图上，默认数组为183.11 MiB，nativeRaster为22.89 MiB。
4. 原组合PNG画布24×22英寸、300 dpi，共4752万像素；单份RGBA像素缓冲约181.3 MiB，设备和渲染可能使用更多缓冲。R图形引擎和grid显示列表还可能保留已经绘制的对象。

## 修复行为

- 不再读取07绘图RDS，旧文件无需删除、转换或重建。
- 从06富集结果重建需要的曲线；按各比较显著条目的FDR排序取 `curve_top_each`，再按出现次数取前三个ID，维持原选择思路。
- 保留 `cfg$contrast_levels` 的比较顺序及九组组合图；缺少结果时显示占位面板。
- 每个面板单独构建grob、绘制、释放；导出设备关闭图形引擎和grid显示列表。
- 从现有Pathview输出目录发现通路图，以nativeRaster逐张读取，保持纵横比；每页最多九张，超过时分页，保留全部图像。
- 复用当前会话中有效的inputs、deg和enrichment对象，避免重复读取；新会话按需读取01/03/06。保留输入、统计参数、富集背景检查，并核对06记录的03文件MD5。
- 临时对象位于局部作用域，执行后释放；不会清除用户已有全局对象。
- PNG默认限制为2400万像素，必要时降低导出dpi；PDF保留原画布尺寸及曲线矢量元素，Pathview图像仍为栅格。
- 输出 `03_analysis/10_objects_and_summary/09_memory_usage.tsv`，记录R堆、RSS和进程峰值RSS（后两项在Linux可用）。RSS峰值包含同一进程之前的操作。

PNG像素预算只限制输出画布，不是整个R进程的内存上限；单张源PNG解码、统计结果及图形对象仍需要内存。

## 应用补丁

将补丁包上传服务器，在项目顶层执行：

```bash
unzip -o Calicivirus_09_memory_fix_20261019.zip
```

也可以仅用新的 `09_combine_export.R` 替换 `analysis_scripts/09_combine_export.R`。

建议重启R且不恢复旧工作区，以释放原07/08及失败09遗留的对象。将R工作目录设为项目顶层，恢复运行06时使用的自定义cfg参数，然后执行：

```r
source("analysis_scripts/00_config.R")
# 如有自定义分析参数，在此恢复
source("analysis_scripts/09_combine_export.R")
```

如果希望在有限内存下仍保持300 dpi，可以进一步缩小组合图画布；如需保持旧版24×22英寸、300 dpi输出，且资源足够，可在运行09前设置：

```r
cfg$combine_max_pixels <- 50000000
```

更保守的PNG预算：

```r
cfg$combine_max_pixels <- 12000000
```

## 验证范围

本地R 4.3.3、ggplot2 3.4.4、patchwork 1.2.0、png 0.1-8：

- 全部24个R文件语法检查通过。
- 将旧07绘图RDS故意替换为不可读文本，模块09仍完成导出，证明修复不依赖该文件。
- 九个比较，每个曲线由三个子图组成，每个子图20000个模拟点；验证PDF和PNG导出。
- 十二张3000×2000像素RGBA通路测试图，验证两页输出、像素预算及图形设备关闭。
- 注入曲线绘图失败，验证FAILED状态记录和设备释放，之后Pathview与汇总仍可执行。
- 验证临时对象不进入全局工作区；查看合成组合图，确认三列布局、标签与栅格纵横比。
- 正常压力测试的进程峰值RSS约428 MiB；这是本地合成图形测试，并非真实转录组数据的内存上限，也不是与旧版同负载的完整性能对照。

测试使用模拟GSEA曲线接口，不重新计算富集。尚未在用户服务器和真实数据上运行新版09；真实峰值取决于结果规模、原始PNG尺寸、R/包版本和现有工作区。

回归测试：

```bash
Rscript analysis_scripts/tests/combine_export_memory.R
```

参考：R官方grid显示列表说明 https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grid.display.list.html ；png作者文档 https://www.rforge.net/png/docs/reference/readPNG.html 。
