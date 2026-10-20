# main 更新后在超算容器执行

工作流：`.github/workflows/scnet-main.yml`。

每次 push 到 main（包括 PR 合入）提交一个新容器，启动脚本拉取该次提交的完整 SHA，在持久项目目录执行 `SCNET_RUN_COMMAND`。
也可在 Actions 页面手动运行，所选分支必须是 main。没有配置 PR 触发器。

## 当前实现边界

- Actions 成功表示超算接受了创建请求，**不表示 R 分析成功**。实例 ID、提交 SHA、远程日志路径写入 Actions Summary。
- 容器创建时 `useStartScript=true`，使用 `header` 范围执行脚本，无需再调用 execute-script，避免重复执行。
- 默认一个容器、CPU 任务。CPU、内存、镜像、区域和挂载必须按账户实际权限填写；示例规格不是资源可用性承诺。
- 通过 `timeoutLimit` 设置有限的容器自动停止时长。**脚本提前完成不会提前停机**；可在超算控制台提前停止。暂未实现远程结束状态回传、完成即停机和结果下载。
- 所有提交均会创建容器。共享项目使用 `flock` 等待前一任务释放锁；等待也会消耗实例资源和超时额度，突发提交很多时应先暂停工作流或增大间隔。
- 取消 GitHub workflow 不会取消已经创建的超算实例。网络超时可能发生在创建成功后；按 Actions 日志中的实例名查重，再决定是否重新运行。

## 1. GitHub 配置

在仓库 Settings → Secrets and variables → Actions 配置。

Repository secrets：

| 名称 | 内容 |
|---|---|
| `SCNET_USER` | 国家超算平台用户名，签名认证使用 |
| `SCNET_ACCESS_KEY` | 下载授权码中的 Accesskey |
| `SCNET_SECRET_KEY` | 对应 Secretkey |

Repository variables：

| 名称 | 内容 |
|---|---|
| `SCNET_CLUSTER_ID` | 有容器权限的计算区域 ID；不能填平台区域 0 |
| `SCNET_CONTAINER_JSON` | 按 `scripts/scnet/container.example.json` 修改后的完整 JSON |
| `SCNET_WORK_ROOT` | 持久化挂载内目录，如 `/mnt/work/transcriptome` |
| `SCNET_RUN_COMMAND` | 容器内执行的分析命令，见下文 |

AK/SK 只存 GitHub secrets，不写进仓库或启动脚本。Actions 中签名与各区域 token 自动遮蔽。
区域凭证每次通过 AK/SK 申请，从授权区域接口取 `enable=true` 的 `aiUrls`；不硬编码区域域名。
兼容返回 URL 已带 `/ai` 的情况。

JSON 里应填写平台提供的 `imagePath`、`version`、`resourceGroup`，以及你的授权 `sourcePath`。
`ramSize` 单位 MB。`targetPath` 是容器路径；`sourcePath` 是超算共享存储路径，两者不要混用。
`SCNET_WORK_ROOT` 必须在 `targetPath` 内。建议使用专用工作目录，避免人工分析与自动更新相互影响。

## 2. 容器与项目准备

容器镜像需要预装 Bash、base64、Git、util-linux 的 flock、R 和全部分析依赖，并能通过 HTTPS 访问 GitHub。
使用包含所需 R/Bioconductor 包的固定镜像。首次在该镜像中验证环境，避免每次提交时重新安装全部依赖。
当前仓库为公开仓库，远端直接 HTTPS 获取；若改成私有仓库，需在容器预置只读 Git 凭证，当前实现不会把 GitHub token 发送到超算。

持久目录结构：

```text
/mnt/work/transcriptome/
  project/              # Git 工作副本及相对路径输入、结果、_targets 缓存
  runs/<run-id>-<attempt>/
    run.log
    commit
    started_at
    finished_at
    exit_code
```

`project/` 不存在时脚本会自动 clone；若已存在，必须是本仓库 Git 工作副本。
已有未提交的受跟踪文件修改时终止，不执行覆盖。未跟踪文件不会被清理，但与仓库文件同名冲突时 checkout 会失败。
分析输入、结果、`_targets/` 应保持在受跟踪源码之外；在适当的 `.gitignore` 或持久工作副本 `.git/info/exclude` 中排除。
为便于保持缓存，工作目录固定，不为每个 SHA 新建分析目录。

在首次正式运行前，将数据准备到 `project/` 所需的相对路径，或将分析配置改为容器内的绝对输入路径。
同一共享存储也可以额外挂载为 `/mnt/data`，让配置指向该路径；不要把 FASTQ、表达矩阵或结果上传到 GitHub。
不同容器中的运行用户必须对共享目录有一致的读写权限，且存储支持跨容器文件锁。

## 3. 配置执行入口

先用下面命令验证通路（作为变量的完整值填写）：

```bash
printf 'container ready\n'; git rev-parse HEAD; Rscript --version
```

合入 targets 工程后，可改为：

```bash
Rscript -e 'stopifnot(file.exists("_targets.R")); targets::tar_make()'
```

targets 缓存位于持久 `project/_targets/`；后续 `tar_make()` 根据依赖变化增量运行。
本次新增的是超算执行层，没有将 `targets-engineering-20261019` 分支自动合入 main。
如使用原模块入口，必须先将对应模块代码合入 main，再设置：

```bash
Rscript analysis_scripts/run_all.R
```

如需多个步骤，建议把入口写成仓库中的脚本，再将变量设为 `bash scripts/run_analysis.sh`。
启动层用 `bash -euo pipefail -c` 执行命令；不要在此配置后台执行 `&`，否则退出码不能代表分析完成。
实际分析并行度须与容器 CPU/RAM 配额匹配；工作流不调整 R 分析参数。

## 4. 检查结果

Actions Summary 给出容器 ID 和本次日志路径。在超算挂载目录查看 `run.log` 和 `exit_code`：

- `0`：启动层、代码更新和分析命令完成。
- 非零：启动或分析失败，检查 `run.log`。
- 缺少 `exit_code`：可能排队、运行中、启动失败、超时或被中止，需结合超算实例状态判断。

实例运行到时间上限时平台应自动停止；自动停止行为、脚本启动行为及网络连通性必须在你的计算区域做一次真实验证。
当前没有超算凭证，已验证范围为本地离线测试，不能声称已经创建或运行超算容器。

## 5. 离线验证

```bash
python3 -m unittest discover -s tests -p 'test_scnet.py' -v
bash -n scripts/scnet/container_entry.sh
```

本地环境中配置上述变量并补充 GitHub 执行上下文后，可运行 `python3 scripts/scnet/submit.py --dry-run`。
dry-run 不登录、不访问 API、不创建容器。测试使用本地 Git fixture，检查旧 SHA 获取、缓存保留和失败退出码。

## 官方接口依据

- [AK/SK 获取凭证](https://www.scnet.cn/ac/openapi/doc/2.0/api/safecertification/get-user-tokens-aksk.html)
- [获取授权区域](https://www.scnet.cn/ac/openapi/doc/2.0/api/safecertification/get-center-info.html)
- [创建容器](https://www.scnet.cn/ac/openapi/doc/2.0/api/container/create.html)
- [已有容器执行脚本](https://www.scnet.cn/ac/openapi/doc/2.0/api/container/execute.html)
- [容器停止](https://www.scnet.cn/ac/openapi/doc/2.0/api/container/stop.html)
- [GitHub Actions 工作流语法](https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-syntax)
