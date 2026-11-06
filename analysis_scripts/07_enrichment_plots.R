if (!exists("cfg")) stop("请先source analysis_scripts/00_config.R")
require_enrichment_packages()
inputs <- read_module("01_inputs")
deg <- read_module("03_deg")
enrichment <- read_module("06_enrichment")
check_module_inputs(inputs)
check_analysis_settings(deg)
check_enrichment_settings(enrichment, inputs)
out_dir <- module_dir("07_enrichment")
combined_dir <- module_dir("09_combined_figures")
plots <- list()
pathview_files <- list()
global_limit <- max(abs(unlist(lapply(deg$comparisons, function(item) item$all$logFC))), na.rm = TRUE)
for (db in c("GO", "KEGG")) {
  database <- if (db == "GO") "GO_BP" else "KEGG"
  plots[[db]] <- list(ORA_dot = list(), ORA_bar = list(), ORA_upset = list(), ORA_network = list(),
                      GSEA_bar = list(), GSEA_network = list(), GSEA_curves = list())
  for (comparison in cfg$contrast_levels) {
    item <- enrichment$results[[comparison]][[db]]
    if (is.null(item)) next
    dir <- file.path(out_dir, database, comparison, "plots")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    for (direction in c("Up", "Down")) {
      key <- paste(comparison, direction, sep = "_")
      res <- subset_enrichment(item$ORA[[direction]])
      if (is.null(res)) {
        for (type in c("ORA_dot", "ORA_bar", "ORA_upset", "ORA_network")) plots[[db]][[type]][key] <- list(NULL)
        next
      }
      plot_calls <- list(
        ORA_dot = function() plot_ora_dot(res, inputs$public_data, showCategory = cfg$top_terms, title = key, out_dir = NULL),
        ORA_bar = function() plot_ora_geom_col(res, inputs$public_data, showCategory = cfg$top_terms, direction = direction, title = key, out_dir = NULL),
        ORA_upset = function() plot_ora_upset(res, inputs$public_data, n = min(8, nrow(result_table(res))), title = key, out_dir = NULL),
        ORA_network = function() plot_ora_cnet(res, inputs$public_data, showCategory = 5, title = key,
                                               logFC = stats::setNames(deg$comparisons[[comparison]]$all$logFC, deg$comparisons[[comparison]]$all$ENTREZID))
      )
      for (type in names(plot_calls)) {
        p <- attempt_plot("07_enrichment_plots", paste(database, type, key), plot_calls[[type]]())
        plots[[db]][[type]][key] <- list(p)
        save_plot_pair(p, file.path(dir, paste0(type, "_", direction)), n_terms = if (type %in% c("ORA_dot", "ORA_bar")) cfg$top_terms else NULL)
      }
    }
    gsea <- subset_enrichment(item$GSEA)
    if (!is.null(gsea)) {
      p <- attempt_plot("07_enrichment_plots", paste(database, comparison, "GSEA_bar"), {
        plot_gsea_double_bar(gsea, inputs$public_data, contrast_name = comparison)
      })
      plots[[db]]$GSEA_bar[comparison] <- list(p)
      save_plot_pair(p, file.path(dir, "GSEA_double_bar"), n_terms = min(cfg$top_terms * 2, nrow(result_table(gsea))))
      p <- attempt_plot("07_enrichment_plots", paste(database, comparison, "GSEA_network"), {
        plot_gsea_cnet(gsea, inputs$public_data, contrast_name = comparison, showCategory = 5, out_dir = NULL)
      })
      plots[[db]]$GSEA_network[comparison] <- list(p)
      save_plot_pair(p, file.path(dir, "GSEA_network"), 12, 9)
      ranked <- result_table(gsea)
      ranked <- ranked[order(ranked$p.adjust), , drop = FALSE]
      ids <- utils::head(ranked$ID, cfg$curve_top_each)
      curves <- attempt_plot("07_enrichment_plots", paste(database, comparison, "GSEA_curves"), {
        plot_gsea_gseaplot2(subset_enrichment(gsea, ids), out_dir = NULL)
      })
      if (!is.null(curves)) {
        plots[[db]]$GSEA_curves[[comparison]] <- curves
        for (id in names(curves)) save_plot_pair(curves[[id]], file.path(dir, paste0("GSEA_curve_", safe_name(id))), 9, 7)
        combine_plots(curves, file.path(dir, "GSEA_representative_curves"), ncol = min(3, length(curves)))
      }
    } else {
      plots[[db]]$GSEA_bar[comparison] <- list(NULL)
      plots[[db]]$GSEA_network[comparison] <- list(NULL)
    }
    if (db == "KEGG") {
      if (cfg$run_pathview && requireNamespace("pathview", quietly = TRUE)) {
        pathview_files[[comparison]] <- list()
        objects <- list(ORA_Up = item$ORA$Up, ORA_Down = item$ORA$Down, GSEA = item$GSEA)
        for (method in names(objects)) {
          files <- attempt_plot("07_pathview", paste(comparison, method), {
            plot_pathview_result(objects[[method]], item$ranking, file.path(dir, "pathview", method), global_limit)
          })
          pathview_files[[comparison]][method] <- list(files)
        }
      } else record_status("07_pathview", comparison, "SKIPPED", "run_pathview=FALSE或缺少pathview包")
    }
  }
  for (time in cfg$time_levels) {
    comparisons <- cfg$contrast_levels[endsWith(cfg$contrast_levels, paste0("NC_", time))]
    keys <- as.vector(t(outer(comparisons, c("Up", "Down"), paste, sep = "_")))
    for (type in c("ORA_dot", "ORA_bar", "ORA_upset", "ORA_network")) {
      panels <- stats::setNames(lapply(keys, function(key) plots[[db]][[type]][[key]]), keys)
      combine_plots(panels, file.path(combined_dir, paste(database, type, time, sep = "_")), 2, 23, 23)
    }
  }
  for (type in c("GSEA_bar", "GSEA_network")) {
    panels <- stats::setNames(lapply(cfg$contrast_levels, function(key) plots[[db]][[type]][[key]]), cfg$contrast_levels)
    combine_plots(panels, file.path(combined_dir, paste(database, type, "nine_contrasts", sep = "_")), 3, 24, 24)
  }
}
save_module(list(plots = plots, pathview_files = pathview_files), "07_enrichment_plots")
