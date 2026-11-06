if (!exists("cfg")) stop("请先source analysis_scripts/00_config.R")

# 局部执行，释放临时对象
local({
  memory_log <- data.frame()
  report_memory <- function(stage) {
    heap <- sum(gc()[, 2])
    rss <- peak <- NA_real_
    if (file.exists("/proc/self/status")) {
      lines <- readLines("/proc/self/status", warn = FALSE)
      value <- function(key) {
        line <- grep(paste0("^", key, ":"), lines, value = TRUE)
        if (length(line)) as.numeric(gsub("[^0-9]", "", line[1])) / 1024 else NA_real_
      }
      rss <- value("VmRSS")
      peak <- value("VmHWM")
    }
    row <- data.frame(Stage = stage, RHeapMiB = heap, RSSMiB = rss, PeakRSSMiB = peak)
    message("09 ", stage, "；R堆=", round(heap), " MiB；RSS=", round(rss), " MiB")
    row
  }
  memory_log <- rbind(memory_log, report_memory("start"))
  inputs <- get0("inputs", inherits = TRUE)
  if (!is.list(inputs) || is.null(inputs$signature)) inputs <- read_module("01_inputs")
  check_module_inputs(inputs)
  deg <- get0("deg", inherits = TRUE)
  if (!is.list(deg) || is.null(deg$settings)) deg <- read_module("03_deg")
  check_analysis_settings(deg)
  deg <- NULL
  object_dir <- module_dir("10_objects_and_summary")
  combined_dir <- module_dir("09_combined_figures")
  enrich_file <- file.path(object_dir, "06_enrichment.rds")
  max_pixels <- if (is.null(cfg$combine_max_pixels)) 24000000 else cfg$combine_max_pixels
  if (length(max_pixels) != 1 || !is.finite(max_pixels) || max_pixels < 1) stop("combine_max_pixels必须为正数")
  if (length(cfg$dpi) != 1 || !is.finite(cfg$dpi) || cfg$dpi < 1) stop("dpi必须为正数")

  # 每次只绘制一个面板
  save_streamed_page <- function(draw_panel, n_panels, prefix, ncol = 3L, width = 24, height = 22) {
    if (!n_panels) return(NULL)
    dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
    nrow <- ceiling(n_panels / ncol)
    png_dpi <- min(cfg$dpi, floor(sqrt(max_pixels / (width * height))))
    if (png_dpi < 1) stop("组合图尺寸超过像素预算，请提高combine_max_pixels")
    if (png_dpi < cfg$dpi) message("PNG像素预算：", basename(prefix), " 使用 ", png_dpi, " dpi；PDF保留原尺寸")
    for (format in c("pdf", "png")) {
      file <- paste0(prefix, ".", format)
      message("09 导出：", basename(file))
      if (format == "pdf") {
        device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
        device(file, width = width, height = height)
      } else {
        args <- list(filename = file, width = width, height = height, units = "in", res = png_dpi)
        if (capabilities("cairo")) args$type <- "cairo"
        do.call(grDevices::png, args)
      }
      device_id <- grDevices::dev.cur()
      tryCatch({
        grDevices::dev.control(displaylist = "inhibit")
        grid::grid.newpage()
        grid::grid.display.list(FALSE)
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))
        for (i in seq_len(n_panels)) {
          grid::pushViewport(grid::viewport(layout.pos.row = (i - 1L) %/% ncol + 1L, layout.pos.col = (i - 1L) %% ncol + 1L))
          tryCatch(draw_panel(i), finally = grid::popViewport())
          invisible(gc())
        }
        grid::popViewport()
      }, finally = grDevices::dev.off(device_id))
    }
    invisible(prefix)
  }

  if (file.exists(enrich_file)) {
    enrichment <- get0("enrichment", inherits = TRUE)
    if (!is.list(enrichment) || is.null(enrichment$signature)) enrichment <- readRDS(enrich_file)
    check_enrichment_settings(enrichment, inputs)
    deg_file <- file.path(object_dir, "03_deg.rds")
    if (!identical(enrichment$signature$deg, unname(tools::md5sum(deg_file)))) stop("差异结果已改变，请重新执行06模块")
    memory_log <- rbind(memory_log, report_memory("results_ready"))
    for (db in c("GO", "KEGG")) {
      # 保留原来的曲线选择规则
      curve_ids <- unlist(lapply(cfg$contrast_levels, function(comparison) {
        tbl <- result_table(enrichment$results[[comparison]][[db]]$GSEA, TRUE)
        if (!nrow(tbl)) return(character())
        utils::head(tbl$ID[order(tbl$p.adjust)], cfg$curve_top_each)
      }), use.names = FALSE)
      ids <- utils::head(names(sort(table(curve_ids), decreasing = TRUE)), 3)
      for (id in ids) {
        draw_curve <- function(i) {
          comparison <- cfg$contrast_levels[i]
          item <- enrichment$results[[comparison]][[db]]
          res <- if (is.null(item)) NULL else subset_enrichment(item$GSEA, id)
          if (is.null(res)) {
            p <- empty_panel_plot(comparison, if (is.null(item)) "Database unavailable" else "Not significant or not tested")
          } else {
            curves <- plot_gsea_gseaplot2(res, pvalue_cutoff = cfg$fdr_cutoff, out_dir = NULL)
            p <- if (length(curves)) as_combination_plot(curves[[1]]) else empty_panel_plot(comparison, "No curve returned")
            title <- paste(comparison, id)
            p <- if (inherits(p, "patchwork")) p + patchwork::plot_annotation(title = title) else p + ggplot2::labs(title = title)
          }
          grob <- if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else ggplot2::ggplotGrob(p)
          grid::grid.draw(grob)
          invisible(NULL)
        }
        attempt_plot("09_combination", paste(db, id, "curves_across_comparisons"), {
          save_streamed_page(draw_curve, length(cfg$contrast_levels),
                             file.path(combined_dir, paste0("14_", db, "_", safe_name(id), "_curves_across_comparisons")))
        })
        memory_log <- rbind(memory_log, report_memory(paste(db, id)))
      }
    }
    enrichment <- NULL
    invisible(gc())
  }

  # 从磁盘发现通路图
  pathview_root <- file.path(cfg$output_dir, "07_enrichment", "KEGG")
  pathview_files <- list.files(pathview_root, pattern = "pathview.*\\.png$", recursive = TRUE, full.names = TRUE)
  pathview_files <- pathview_files[grepl("/pathview/", gsub("\\\\", "/", pathview_files), fixed = TRUE)]
  if (length(pathview_files) && requireNamespace("png", quietly = TRUE)) {
    labels <- substring(dirname(pathview_files), nchar(pathview_root) + 2L)
    comparisons <- sub("/.*", "", labels)
    file_order <- order(match(comparisons, cfg$contrast_levels), labels, basename(pathview_files))
    pathview_files <- pathview_files[file_order]
    for (name in unique(basename(pathview_files))) {
      files <- pathview_files[basename(pathview_files) == name]
      # 每页最多九张，保留所有通路图
      pages <- split(seq_along(files), ceiling(seq_along(files) / 9L))
      for (page in seq_along(pages)) {
        page_files <- files[pages[[page]]]
        draw_pathview <- function(i) {
          file <- page_files[i]
          image <- png::readPNG(file, native = TRUE)
          label <- substring(dirname(file), nchar(pathview_root) + 2L)
          grid::grid.text(label, x = 0.5, y = 0.98, just = "top", gp = grid::gpar(fontsize = 10))
          ratio <- ncol(image) / nrow(image)
          panel_ratio <- grid::convertWidth(grid::unit(1, "npc"), "in", valueOnly = TRUE) /
            grid::convertHeight(grid::unit(0.9, "npc"), "in", valueOnly = TRUE)
          w <- min(1, ratio / panel_ratio)
          h <- 0.9 * min(1, panel_ratio / ratio)
          grid::grid.raster(image, x = 0.5, y = 0.46, width = w, height = h, interpolate = FALSE)
          invisible(NULL)
        }
        suffix <- if (length(pages) > 1) sprintf("_page%02d", page) else ""
        prefix <- file.path(combined_dir, paste0("15_pathview_RASTER_", safe_name(name), suffix))
        attempt_plot("09_combination", basename(prefix), {
          save_streamed_page(draw_pathview, length(page_files), prefix, height = ceiling(length(page_files) / 3) * 7)
        })
        memory_log <- rbind(memory_log, report_memory(paste("Pathview", name, page)))
      }
    }
  } else if (length(pathview_files)) record_status("09_combination", "Pathview", "SKIPPED", "缺少png包")

  # 汇总现有输出及状态
  files <- list.files(cfg$output_dir, recursive = TRUE, full.names = TRUE)
  figure_files <- files[grepl("\\.(pdf|png)$", files, ignore.case = TRUE)]
  inventory <- data.frame(Path = sub(paste0(cfg$output_dir, "/"), "", figure_files, fixed = TRUE),
                          Format = tolower(tools::file_ext(figure_files)), SizeBytes = file.info(figure_files)$size)
  save_table(inventory, file.path(object_dir, "figure_inventory.tsv"))
  capture.output(sessionInfo(), file = file.path(object_dir, "sessionInfo.txt"))
  settings <- data.frame(Parameter = names(cfg), Value = vapply(cfg, function(x) paste(x, collapse = ";"), character(1)))
  save_table(settings, file.path(object_dir, "analysis_parameters.tsv"))
  record_status("09_summary", "output_inventory", "SUCCESS", paste(nrow(inventory), "figure files"))
  statuses <- data.table::fread(file.path(object_dir, "module_status.tsv"), quote = "", data.table = FALSE)
  summary <- c("# 分析运行汇总", "", paste0("样本数：", nrow(inputs$sample_info)), paste0("比较数：", nrow(inputs$contrasts)),
               paste0("差异方法：", cfg$method), paste0("显著性：FDR < ", cfg$fdr_cutoff, "；|log2FC| >= ", cfg$abs_logfc_cutoff),
               paste0("参考注释：", cfg$reference_species, "；tax_id=", cfg$tax_id), paste0("GO：BP；祖先传播=", cfg$go_ancestors),
               paste0("图像文件数：", nrow(inventory)), "", "## 状态计数", "", capture.output(print(table(statuses$Status))), "",
               "请查看 module_status.tsv 中的 FAILED、SKIPPED、NO_SIGNIFICANT_RESULT、NOT_CONFIGURED；空结果不代表代码失败。",
               "Pathview组合图包含栅格通路图。ORA方向表示显著DEG集合；GSEA方向表示NES符号。",
               paste0("09逐面板导出；PNG像素预算=", max_pixels, "，PDF保留原尺寸；Pathview每页最多九张。"))
  writeLines(summary, file.path(object_dir, "run_summary.md"))
  memory_log <- rbind(memory_log, report_memory("complete"))
  save_table(memory_log, file.path(object_dir, "09_memory_usage.tsv"))
  message("汇总完成：", cfg$output_dir)
  if (any(statuses$Status == "FAILED")) warning("部分背景或图像失败，请查看module_status.tsv；已完成结果已保存", call. = FALSE)
})
