console_app <- function(code_dir, data_dir = code_dir, config_dir = "config") {
  context <- console_context_lst(code_dir, data_dir, config_dir)
  initial <- console_read_config_lst(context$config_dir)
  controller <- console_new_controller()
  available_cores <- console_available_cores()
  has_graph <- requireNamespace("visNetwork", quietly = TRUE)
  modules <- console_modules_lst()
  module_choices <- stats::setNames(names(modules), vapply(modules, `[[`, character(1), "label"))
  ui <- shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(shiny::HTML(paste0(
      "body{background:#f5f7fa;color:#18263a;font-family:system-ui,sans-serif}",
      ".container-fluid{max-width:1500px;margin:auto;padding:24px 30px}",
      ".hero{display:flex;justify-content:space-between;align-items:center;margin:8px 0 25px}",
      ".hero h2{font-weight:700;margin:0 0 8px}.muted{color:#63738a}.label{letter-spacing:.6px}",
      ".well,.tab-content{background:white;border:1px solid #dde4ec;border-radius:10px;box-shadow:none}",
      ".tab-content{padding:24px;margin-top:16px}.well{padding:20px}.btn-primary{background:#21635d;border-color:#21635d}",
      "pre{background:#122231;color:#d8e8e9;border:0;max-height:480px;overflow:auto;white-space:pre-wrap}",
      "table{font-size:13px}.result-image img{max-width:100%;height:auto}.path{word-break:break-all}",
      ".nav-tabs>li>a{font-size:15px}.btn{margin:3px 4px 3px 0}.status-box{padding:14px;background:#edf4f4;border-radius:8px;margin-bottom:18px}"
    ))), shiny::tags$script(shiny::HTML(
      "Shiny.addCustomMessageHandler('consoleBusy',function(x){document.querySelectorAll('.freeze input,.freeze select,.freeze button').forEach(function(e){e.disabled=x;});});"
    ))),
    shiny::tags$div(class = "hero", shiny::tags$div(shiny::tags$h2("转录组分析控制台"),
      shiny::tags$div(class = "muted", "RStudio Server · targets · 试验验证版")), shiny::tags$span(class = "label label-default", "MVP")),
    shiny::uiOutput("status"),
    shiny::tabsetPanel(id = "page",
      shiny::tabPanel("1 · 项目与输入", value = "inputs",
        shiny::tags$p(class = "muted path", paste("代码目录：", context$code_dir)),
        shiny::tags$div(class = "freeze",
          shiny::fluidRow(shiny::column(6, shiny::textInput("data_dir", "服务器数据目录", context$data_dir)),
            shiny::column(6, shiny::textInput("config_dir", "配置目录", context$config_dir))),
          shiny::actionButton("load_config", "读取配置"),
          shiny::fluidRow(shiny::column(6, shiny::textInput("count_file", "Count 表", initial$paths$count_file),
            shiny::textInput("gtf_file", "GTF 注释", initial$paths$gtf_file)),
            shiny::column(6, shiny::textInput("sample_file", "样本元数据表", initial$paths$sample_file),
              shiny::textInput("gene2go_file", "gene2go 注释", initial$paths$gene2go_file))),
          shiny::textInput("output_dir", "结果目录（相对路径基于数据目录）", initial$paths$output_dir),
          shiny::fluidRow(shiny::column(4, shiny::numericInput("fdr", "FDR 阈值", initial$analysis$filter$fdr, min = 0.000001, max = 1, step = 0.01)),
            shiny::column(4, shiny::numericInput("abs_logfc", "|log2FC| 阈值", initial$analysis$filter$abs_logfc, min = 0, step = 0.1)),
            shiny::column(4, shiny::numericInput("top_terms", "富集单图 Top N", initial$figures$display$top_terms, min = 1, step = 1))),
          shiny::actionButton("save_config", "保存配置快照"), shiny::actionButton("check", "检查输入与环境", class = "btn-primary")
        ),
        shiny::tags$p(class = "muted", "界面修改保存为独立快照。分组、比较及其他参数沿用所选配置目录；首版通过 YAML 编辑。"),
        shiny::verbatimTextOutput("config_note"), shiny::tableOutput("input_table")),
      shiny::tabPanel("2 · 样本与比较", value = "design",
        shiny::tags$h4("样本元数据"), shiny::tableOutput("metadata"),
        shiny::tags$h4("将执行的比较"), shiny::tableOutput("comparisons"),
        shiny::tags$h4("环境依赖概览"), shiny::tableOutput("packages"),
        shiny::tags$p(class = "muted", "表格各显示前100行。输入检查只读取count表头；分析阶段继续验证完整计数和注释。完整依赖见 scripts/check_environment.R。")),
      shiny::tabPanel("3 · 运行与日志", value = "run",
        shiny::fluidRow(shiny::column(4, shiny::wellPanel(shiny::tags$div(class = "freeze",
          shiny::checkboxInput("run_all", "运行完整工程（包含汇总）", FALSE),
          shiny::checkboxGroupInput("modules", "选择模块", choices = module_choices, selected = "qc"),
          shiny::checkboxInput("manual_selection", "高级：按target名称选择", FALSE),
          shiny::selectizeInput("target_names", "指定targets（先生成依赖预览）", choices = character(), multiple = TRUE),
          shiny::tags$h4("targets 执行选项"),
          shiny::tags$p(class = "muted", paste("检测到可用CPU：", available_cores, "；worker数×线程数不得超过此值。")),
          shiny::numericInput("workers", "并行worker数（独立R进程）", 1, min = 1, max = available_cores, step = 1),
          shiny::numericInput("threads_per_worker", "每worker线程上限", 1, min = 1, max = available_cores, step = 1),
          shiny::tags$p(class = "muted", "并行worker使用crew。线程上限作用于data.table及支持环境变量的BLAS/OpenMP库；GSEA内部仍串行。"),
          shiny::selectInput("reporter", "日志级别", c("逐任务" = "verbose", "带时间戳" = "timestamp", "平衡" = "balanced",
            "精简" = "terse", "静默" = "silent")),
          shiny::selectInput("error", "错误策略", c("遇错停止" = "stop", "继续独立任务" = "continue")),
          shiny::selectInput("memory", "内存策略", c("及时释放" = "transient", "自动" = "auto", "保留对象" = "persistent")),
          shiny::numericInput("garbage_collection", "每多少个target执行一次GC（0关闭）", 0, min = 0, step = 1),
          shiny::checkboxInput("shortcut", "shortcut：假定上游已更新，跳过上游检查", FALSE),
          shiny::actionButton("start", "启动分析", class = "btn-primary")),
          shiny::actionButton("stop", "停止当前任务", class = "btn-danger"),
          shiny::tags$p(class = "muted", "模块运行包含必要的上游依赖。停止会保留已完成缓存；再次启动时检查并补算。"))),
          shiny::column(8, shiny::tags$h4("实时日志 · 最近200行"), shiny::verbatimTextOutput("log"),
            shiny::tags$h4("crew worker 日志"), shiny::selectInput("crew_file", "选择worker日志", character()),
            shiny::verbatimTextOutput("crew_log"))),
        shiny::tags$h4("targets 任务状态"), shiny::tableOutput("progress"),
        shiny::tags$p(class = "muted", "状态表最多显示200项，不将任务计数换算为剩余运行时间。")),
      shiny::tabPanel("4 · 结果浏览", value = "results",
        shiny::tags$p(class = "muted", "按需读取单个结果；PNG预览上限20 MiB，表格预览100行，目录最多列出3000个文件。"),
        shiny::actionButton("refresh_results", "刷新结果清单"),
        shiny::fluidRow(shiny::column(4, shiny::selectInput("format", "文件类型", c("全部", "png", "pdf", "tsv", "csv", "txt", "md")),
          shiny::textInput("search", "按路径筛选"), shiny::selectInput("result_file", "选择文件", character()),
          shiny::downloadButton("download", "下载所选文件")),
          shiny::column(8, shiny::textOutput("result_path"), shiny::uiOutput("preview_ui")))),
      shiny::tabPanel("5 · 依赖关系与运行预览", value = "graph",
        shiny::tags$div(class = "freeze", shiny::checkboxInput("graph_all", "展示整个工程（否则使用运行页的目标选择）", TRUE),
          shiny::checkboxInput("graph_outdated", "检查待更新targets（可能读取文件校验信息）", FALSE),
          shiny::actionButton("refresh_graph", "生成依赖预览", class = "btn-primary")),
        shiny::tags$p(class = "muted", "预览在独立进程生成，不执行分析target。图像为生成时的状态快照，可缩放并点击节点查看关联。"),
        shiny::textOutput("plan_note"),
        if (has_graph) visNetwork::visNetworkOutput("dependency_graph", height = "700px") else shiny::tags$p("请安装visNetwork后重启控制台。"),
        shiny::tags$h4("target清单（前200项）"), shiny::tableOutput("manifest"))
    )
  )
  server <- function(input, output, session) {
    app_controller <- controller
    plan_value <- shiny::reactiveVal(NULL)
    crew_choices <- shiny::reactiveVal(character())
    config_note <- shiny::reactiveVal("尚未保存配置或检查输入。")
    result_root <- shiny::reactiveVal(NULL)
    result_files <- shiny::reactiveVal(data.frame(Path = character(), Format = character(), Bytes = numeric()))
    tick <- shiny::reactive({ shiny::invalidateLater(1000, session); console_poll_lst(controller) })
    current_request <- function() {
      ctx <- console_context_lst(context$code_dir, input$data_dir, input$config_dir)
      fields <- c("count_file", "sample_file", "gtf_file", "gene2go_file", "output_dir", "fdr", "abs_logfc", "top_terms")
      overrides <- stats::setNames(lapply(fields, function(field) input[[field]]), fields)
      execution_fields <- c("workers", "threads_per_worker", "reporter", "shortcut", "memory", "garbage_collection", "error", "graph_outdated")
      execution <- stats::setNames(lapply(execution_fields, function(field) input[[field]]), execution_fields)
      execution <- Filter(Negate(is.null), execution)
      list(context = ctx, config = console_effective_config_lst(ctx, overrides), execution = console_execution_lst(execution))
    }
    notify_error <- function(e) shiny::showNotification(conditionMessage(e), type = "error", duration = 10)
    shiny::observe({
      status <- tick()
      session$sendCustomMessage("consoleBusy", identical(status$state, "RUNNING"))
      if (!is.null(controller$report)) result_root(controller$report$output_dir)
      if (!is.null(controller$plan)) plan_value(controller$plan)
      files <- if (is.null(controller$run_dir)) character() else {
        head(list.files(file.path(controller$run_dir, "crew"), recursive = TRUE), 100L)
      }
      if (!identical(files, shiny::isolate(crew_choices()))) {
        crew_choices(files)
        selected <- shiny::isolate(input$crew_file)
        if (!length(selected) || !selected %in% files) selected <- head(files, 1L)
        shiny::updateSelectInput(session, "crew_file", choices = files, selected = selected)
      }
    })
    output$status <- shiny::renderUI({
      status <- tick()
      labels <- c(IDLE = "就绪", RUNNING = "运行中", SUCCEEDED = "已完成", FAILED = "失败", STOPPED = "已停止")
      shiny::tags$div(class = "status-box", shiny::tags$strong(unname(labels[status$state])),
        if (!is.null(status$run_dir)) shiny::tags$span(class = "path", paste0(" · PID ", status$pid, " · ", basename(status$run_dir))),
        if (nzchar(status$error)) shiny::tags$p(status$error))
    })
    shiny::observeEvent(input$load_config, {
      tryCatch({
        if (tick()$state == "RUNNING") stop("请在任务结束后读取另一套配置")
        ctx <- console_context_lst(context$code_dir, input$data_dir, input$config_dir)
        config <- console_read_config_lst(ctx$config_dir)
        for (field in names(config$paths)) shiny::updateTextInput(session, field, value = config$paths[[field]])
        shiny::updateNumericInput(session, "fdr", value = config$analysis$filter$fdr)
        shiny::updateNumericInput(session, "abs_logfc", value = config$analysis$filter$abs_logfc)
        shiny::updateNumericInput(session, "top_terms", value = config$figures$display$top_terms)
        config_note("已读取配置；可编辑路径和阈值，再检查输入。")
      }, error = notify_error)
    })
    shiny::observeEvent(input$save_config, {
      tryCatch({
        if (tick()$state == "RUNNING") stop("运行中不能替换配置快照")
        request <- current_request()
        config_note(paste("配置快照：", console_snapshot_path(request$context, request$config)))
      }, error = notify_error)
    })
    start <- function(action) {
      tryCatch({
        request <- current_request()
        all <- isTRUE(input$run_all) || (action == "plan" && isTRUE(input$graph_all))
        selected <- if (!all && isTRUE(input$manual_selection)) as.character(input$target_names) else NULL
        run_dir <- console_start_run(controller, request$context, request$config, action, input$modules, all,
          execution = request$execution, target_names = selected)
        config_note(paste("本次记录：", run_dir))
        shiny::updateTabsetPanel(session, "page", selected = if (action == "plan") "graph" else "run")
      }, error = notify_error)
    }
    shiny::observeEvent(input$check, start("check"))
    shiny::observeEvent(input$start, start("build"))
    shiny::observeEvent(input$refresh_graph, start("plan"))
    shiny::observeEvent(input$stop, tryCatch(console_stop_run(controller), error = notify_error))
    output$config_note <- shiny::renderText(config_note())
    report <- shiny::reactive({ tick(); controller$report })
    output$input_table <- shiny::renderTable({ shiny::req(report()); report()$inputs })
    output$metadata <- shiny::renderTable({ shiny::req(report()); head(report()$metadata, 100L) })
    output$comparisons <- shiny::renderTable({ shiny::req(report()); head(report()$comparisons, 100L) })
    output$packages <- shiny::renderTable({ shiny::req(report()); report()$packages })
    output$log <- shiny::renderText({
      tick()
      if (is.null(controller$run_dir)) "选择模块后启动分析，或先检查输入。" else console_log_tail(file.path(controller$run_dir, "run.log"))
    })
    output$crew_log <- shiny::renderText({
      tick()
      if (is.null(controller$run_dir) || !length(input$crew_file) || !nzchar(input$crew_file)) return("暂无worker日志；串行运行不启动crew。")
      root <- file.path(controller$run_dir, "crew")
      if (!dir.exists(root)) return("暂无worker日志。")
      console_log_tail(console_result_path(root, input$crew_file))
    })
    shiny::observeEvent(plan_value(), {
      choices <- plan_value()$manifest$name
      selected <- intersect(shiny::isolate(input$target_names), choices)
      shiny::updateSelectizeInput(session, "target_names", choices = choices, selected = selected, server = TRUE)
    })
    output$plan_note <- shiny::renderText({
      plan <- plan_value()
      if (is.null(plan)) return("尚未生成依赖预览。")
      paste("生成于", format(plan$created), "；图中target数：", sum(plan$manifest$In_graph),
        "；待更新数：", if (is.null(plan$outdated)) "未检查" else length(plan$outdated))
    })
    output$manifest <- shiny::renderTable({
      shiny::req(plan_value())
      head(plan_value()$manifest[, c("name", "Selected", "In_graph", "Outdated"), drop = FALSE], 200L)
    })
    if (has_graph) output$dependency_graph <- visNetwork::renderVisNetwork({ shiny::req(plan_value()); plan_value()$graph })
    output$progress <- shiny::renderTable({
      tick()
      if (!requireNamespace("targets", quietly = TRUE) || !dir.exists(context$store)) return(NULL)
      tryCatch(head(as.data.frame(targets::tar_progress(store = context$store)), 200L), error = function(e) NULL)
    })
    shiny::observeEvent(input$refresh_results, {
      tryCatch({
        if (is.null(controller$report)) {
          request <- current_request()
          result_root(console_absolute_path(request$config$paths$output_dir, request$context$data_dir, must_work = FALSE))
        } else result_root(controller$report$output_dir)
        result_files(console_result_files_tbl(result_root()))
      }, error = notify_error)
    })
    shiny::observe({
      files <- result_files()
      if (!is.null(input$format) && input$format != "全部") files <- files[files$Format == input$format, , drop = FALSE]
      if (!is.null(input$search) && nzchar(input$search)) files <- files[grepl(input$search, files$Path, fixed = TRUE), , drop = FALSE]
      selected <- shiny::isolate(input$result_file)
      if (!length(selected) || !selected %in% files$Path) selected <- head(files$Path, 1L)
      shiny::updateSelectInput(session, "result_file", choices = files$Path, selected = selected)
    })
    selected_file <- shiny::reactive({
      shiny::req(result_root(), input$result_file)
      console_result_path(result_root(), input$result_file)
    })
    output$result_path <- shiny::renderText(selected_file())
    output$preview_ui <- shiny::renderUI({
      file <- selected_file()
      ext <- tolower(tools::file_ext(file))
      if (ext == "png" && file.info(file)$size <= 20 * 1024^2) return(shiny::tags$div(class = "result-image", shiny::imageOutput("image", height = "auto")))
      if (ext %in% c("tsv", "csv")) return(shiny::tableOutput("table"))
      if (ext %in% c("txt", "md")) return(shiny::verbatimTextOutput("text"))
      shiny::tags$p("请下载查看该文件。")
    })
    output$image <- shiny::renderImage({
      file <- selected_file()
      shiny::req(tolower(tools::file_ext(file)) == "png", file.info(file)$size <= 20 * 1024^2)
      list(src = file, contentType = "image/png", alt = basename(file))
    }, deleteFile = FALSE)
    output$table <- shiny::renderTable({
      file <- selected_file()
      shiny::req(tolower(tools::file_ext(file)) %in% c("tsv", "csv"))
      head(data.table::fread(file, nrows = 100L, data.table = FALSE), 100L)
    })
    output$text <- shiny::renderText({
      file <- selected_file()
      shiny::req(tolower(tools::file_ext(file)) %in% c("txt", "md"))
      console_log_tail(file)
    })
    output$download <- shiny::downloadHandler(filename = function() basename(selected_file()), content = function(file) {
      if (!file.copy(selected_file(), file, overwrite = TRUE)) stop("下载文件复制失败")
    })
  }
  shiny::onStop(function() console_stop_run(controller))
  shiny::shinyApp(ui, server)
}
