# RStudio Server 控制台 · 试验版

适用于远程Linux服务器上的RStudio Server。默认界面和分析分别运行在独立R进程中，RStudio Console可继续使用；前台对照模式在主会话运行界面。界面仅监听服务器127.0.0.1，由RStudio代理打开，不需要开放新的公网端口。计算使用当前服务器/容器中的资源。

## 首次启动

在代码根目录、已安装原工程分析依赖的R会话中执行：

```r
install.packages(c("shiny", "callr", "yaml", "filelock", "digest", "data.table", "tibble", "dplyr", "rstudioapi"))
source("scripts/console.R")
h <- launch_console(port = 42418)
```

默认使用RStudio Viewer，传入本地地址，由IDE处理端口代理；可从Viewer弹出页面。`viewer = FALSE`通过`rstudioapi::translateLocalUrl()`生成完整地址并尝试打开浏览器。不要直接在自己电脑的浏览器打开服务器localhost地址。启动信息只确认本地端口已监听，不能证明超算网关已成功转发到浏览器。

## 更新代码后重新激活

在包含`_targets.R`的代码根目录更新试验分支，然后在同一个RStudio会话执行：

```r
# 关闭旧面板及其任务
if (exists("stop_console", mode = "function")) stop_console()
source("scripts/console.R")
stopifnot("port" %in% names(formals(launch_console)))
h <- launch_console(port = 42418, viewer = TRUE)
```

`git pull`或覆盖文件不会自动更新当前R内存中的旧函数，必须重新`source()`。若仍提示`unused argument (port = ...)`，检查`getwd()`是否为正确代码目录，以及该文件是否已经更新；`names(formals(launch_console))`应包含port和background。

| 参数 | 默认值 | 含义 |
|---|---|---|
| port | `getOption("transcriptome.console.port", NULL)` | 1–65535整数；NULL随机选择端口；非root用户通常应使用大于1023的端口 |
| background | TRUE | TRUE独立进程运行界面；FALSE在当前R会话以前台方式运行界面 |
| viewer | TRUE | RStudio中使用Viewer；FALSE尝试打开完整代理URL |
| open | TRUE | FALSE只启动服务并返回句柄，不打开页面 |

端口被其他服务占用时直接报错。再次启动同一代码工程且端口相同，会复用现有进程；更换端口或切换到前台模式需显式`stop_console()`。重新加载入口也可复用旧版句柄。

## 网关故障时的前台对照

如果后台本地HTTP访问成功，但Viewer仍出现`Host not found`等错误，先停止后台控制台，再让Shiny直接在RStudio主会话启动相同页面：

```r
stop_console()
launch_console(port = 42418, background = FALSE, viewer = TRUE)
```

这时R Console会被界面占用，使用RStudio Stop按钮结束；界面内的分析任务仍使用独立计算进程。`stop_console()`只管理后台控制台，不用于停止正在运行的前台Shiny循环。固定端口和前台模式不保证解决平台网关问题，但可以避免随机端口并区分后台启动方式与通用转发故障。

后台模式可查看诊断信息：

```r
h <- launch_console(port = 42418, open = FALSE)
h$process$is_alive()
h$port
h$url
h$proxy_path
h$proxy_url
cat(tail(readLines(h$log_file, warn = FALSE), 20), sep = "\n")
```

句柄中的port为实际监听端口（integer），url为本地URL（character），proxy_path/proxy_url为RStudio提供的相对/完整代理地址（character；不在RStudio或转换失败时NULL），process为callr进程对象，log_file为启动日志路径（character）。这些值不能证明外部网关可达。排查时可在RStudio Terminal执行`curl --noproxy '*' --max-time 5 -sS -o /dev/null -w 'HTTP %{http_code}\n' http://127.0.0.1:42418/`。本地HTTP 200仅证明当前环境内可访问，仍需在浏览器验收代理。

代码与数据分开放置时：

```r
source("scripts/console.R")
launch_console(
  code_dir = getwd(),
  data_dir = "/work/home/acx2ojxu3z/zz_杯状病毒转录组数据分析_2026年10月19日",
  config_dir = "config"
)
```

入口也沿用`TRANSCRIPTOME_PROJECT_DIR`和`TRANSCRIPTOME_CONFIG_DIR`环境变量。data_dir/config_dir的相对路径基于代码根目录；count、样本、注释、输出的相对路径基于数据目录。界面不会自动安装包，也不复制主会话的全局分析对象。子进程继承启动时的`.libPaths()`；使用renv时，先在主会话`renv::load()`。

## Addins 菜单入口

在RStudio Terminal、代码根目录执行：

```bash
R CMD INSTALL tools/transcriptomeconsole
```

若默认R包目录不可写，在Console安装到已有的可写个人库：

```r
install.packages("tools/transcriptomeconsole", repos = NULL, type = "source", lib = "/你的可写R包目录")
```

安装后在Addins菜单选择`Transcriptome Console (MVP)`。它从当前目录向上寻找代码根目录；也可设置`TRANSCRIPTOME_CODE_DIR`。菜单`Stop Transcriptome Console`关闭面板及后台任务。Addin包只是入口，实际控制台代码仍由试验分支维护；更新代码后关闭并重新打开面板。若希望Addin固定使用42418端口，在当前R会话先执行`options(transcriptome.console.port = 42418)`；默认NULL仍随机选择。

## 按页面操作

1. **项目与输入**：填写服务器上的数据/配置目录。更换配置目录后点击“读取配置”；编辑五个路径、FDR、|log2FC|和富集单图Top N。点击“保存配置快照”可保存当前设置；启动任务时也自动保存。
2. **检查输入与环境**：后台读取count表头，使用原工程元数据与比较函数检查样本对应、变量类型、比较组对和重复数。完成后在“样本与比较”查看结果。完整count数值、GTF/gene2go内容仍由分析模块验证。环境表只展示依赖概览，完整依赖检查仍使用`scripts/check_environment.R`。
3. **运行与日志**：选择模块并启动；也可勾选“运行完整工程”。运行前再次检查输入。模块运行会检查其上游依赖；例如QC当前依赖注释背景，可能获取KEGG快照，具体沿用原依赖图。“完整工程”包括最终汇总；局部运行不会额外触发整个汇总流程。
4. **停止当前任务**：终止分析进程及其后代，保留已经写入的targets缓存。未完成目标下次重新构建，不提供在函数中间暂停后继续的能力。
5. **结果浏览**：首次或运行后点击“刷新结果清单”，按类型和路径文字筛选，预览单张PNG/表格/文本并下载文件。PDF通过下载查看。清单基于最近一次输入检查的输出目录；若未检查，则使用当前输入框的输出目录。

首版不提供分组/比较表格编辑器；修改`variables.yml`、`comparisons.yml`后，重新检查。原有火山图、MA、rank及QC计算/绘图函数保持原实现。

## 模块与目标

| 模块 | targets目标 |
|---|---|
| 输入与注释 | input_files、annotation_files |
| 完整QC | qc_files |
| 差异结果表 | deg_files |
| 火山图 / MA / rank | deg_plot_files、deg_combined_files |
| 基因集合与表达热图 | gene_set_files、expression_heatmap_files |
| ORA / GSEA结果 | ora_up_files、ora_down_files、gsea_files |
| 富集单图 | ora_plot_files、gsea_plot_files |
| 比较富集与组合 | comparative_files、comparative_combined_files、curve_combined_files |
| Pathview | pathview_files，须在配置中开启 |
| 完整工程 | `_targets.R`的全部目标，包括summary_files |

## 配置、缓存和运行记录

界面不改写仓库中的`config/*.yml`，编辑结果存档到`.console/configs/<SHA256>/`，相同内容复用存档。启动任务前，在持有工程锁时将配置准备到固定的`.console/active-config/`，只替换内容改变的文件，避免修改绘图参数时因整套配置路径改变而重建其他模块。首次从原config切换到执行配置时，文件路径改变，会触发相关依赖检查/重建；后续可按文件内容增量检查。存档快照按只读使用；若被外部修改，控制台会拒绝复用。

仍使用代码根目录`_targets/`缓存和原有结果目录，不新增另一套分析缓存。运行中固定本次配置；手动修改源YAML不会改变已启动任务的快照。源码和输入数据文件不复制，运行中不要编辑它们。

同一代码工程通过`.console/build.lock`防止多个控制台同时构建；锁由界面进程持有，在后台任务退出/停止时释放。该锁不能阻止用户在其他终端直接运行`tar_make()`，因此不要同时从其他入口构建相同工程。共享文件系统需要支持文件锁；多个代码克隆共用同一输出目录的情况不在首版互斥范围内。

| 路径 | 内容与用途 |
|---|---|
| `.console/app_*.log` | Shiny启动/运行日志，排查页面无法打开 |
| `.console/configs/<SHA256>/*.yml` | 七个有效配置，字段与原config文件一致 |
| `.console/active-config/*.yml` | 固定路径的执行配置，仅在任务启动前更新；运行期间不应手动修改 |
| `.console/runs/<时间>_<随机标识>/request.rds` | 本次操作、目标、代码/数据/配置目录、快照路径及创建时间 |
| 同目录`run.log` | 分析标准输出和错误输出 |
| 同目录`validation.rds` | 输入检查通过后保存的元数据、比较、依赖概览、输入路径和输出目录 |
| 同目录`result.rds` | 状态、错误信息和结束时间，由控制台观察任务结束或停止后写入 |

request.rds为list：`context`为代码/数据/配置/状态/缓存路径，`config_dir`为固定执行配置路径，`snapshot_dir`为本次存档快照路径，`action`为check/build，`targets`为目标向量（NULL表示完整工程），`created`为POSIXct时间。后续运行会更新执行配置，复核历史参数请读取snapshot_dir。

validation.rds为list：`metadata`和`comparisons`沿用原工程字段，见README及变量配置；`packages`每行一个包，`Package`为包名（character），`Installed`为是否找到安装目录（logical，不保证包可加载）；`inputs`每行一个输入，`Input`为配置键（character），`Path`为绝对路径（character）；`output_dir`为输出绝对路径（character），`note`为检查范围说明（character）。

result.rds为list：`state`为SUCCEEDED/FAILED/STOPPED（character），`error`为错误或停止说明（character，成功为空字符串），`finished`为结束时间（POSIXct）。界面另有IDLE/RUNNING状态。进程被强行关闭时可能没有result.rds，应结合run.log判断；运行历史首版从文件保留，不提供界面中的历史任务恢复列表。

结果目录清单在内存中生成，每行一个文件：`Path`为输出目录内相对路径（character），`Format`为小写扩展名（character），`Bytes`为文件字节数（numeric）。界面新增功能不生成TSV/XLSX；分析模块原有输出与字段定义保持原文档。

## 页面与任务生命周期

关闭浏览器页面不会停止独立Shiny进程或正在运行的任务；同一RStudio会话再次`launch_console()`会复用已启动进程。明确关闭：

```r
stop_console()
```

关闭控制台会终止其正在运行的分析。主RStudio R会话退出时，受监督的子进程会被结束；该MVP不承诺任务脱离会话长期运行，也不自动提交Slurm。若只想停止分析但保留面板，使用页面中的“停止当前任务”。

日志每次最多读取64 KiB/200行；只按需预览一张PNG（上限20 MiB），表格最多100行，结果清单最多3000项，targets状态最多200项。结果预览拒绝通过路径或符号链接访问输出目录之外的文件。目录扫描仍遍历结果目录，大量结果建议按路径筛选并适当整理。

## 试验验收

自动验证：

```bash
Rscript tests/run_console_tests.R
```

覆盖配置稳定性、输入/比较检查、模块映射、日志边界、结果目录限制、后台检查、互斥、失败恢复、真实targets最小构建及缓存复用、停止及后代进程清理、Shiny服务端和HTTP启动、固定端口、端口冲突、旧句柄复用、Addin端口选项、前台HTTP对照、Addin安装。最小targets用合成项目验证调度，不替代真实转录组计算或RStudio Server代理验收。测试不需要原始数据或完整生信包；依赖由`.github/workflows/console-mvp.yml`安装。

服务器首次验收依次做：打开面板并确认Console可输入；检查样本及九个默认比较；仅运行QC；检查日志和PNG；再次运行观察缓存跳过；启动一个需重建的模块并停止；重新启动补算。最后验证Addins入口和Viewer备用方式。RStudio Server的实际反向代理、真实数据规模及服务器依赖仍需在你的环境验收。
