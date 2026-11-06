theme_project <- function(style) {
  ggplot2::theme_bw(base_size = style$base_size, base_family = style$base_family) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = style$border_width),
      plot.title = ggplot2::element_text(hjust = 0.5, size = style$title_size), strip.background = ggplot2::element_rect(fill = NA),
      legend.key = ggplot2::element_blank(), axis.text = ggplot2::element_text(color = "black"))
}

theme_pca <- function(style) theme_project(style) + ggplot2::theme(aspect.ratio = 1)

theme_enrichment <- function(style) {
  theme_project(style) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = style$axis_angle, hjust = 1, vjust = 1),
                                      axis.title = ggplot2::element_blank())
}

empty_result_plot <- function(title, detail = "No significant result") {
  ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = detail, size = 4) +
    ggplot2::labs(title = title) + ggplot2::theme_void() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

as_panel_grob <- function(plot) {
  if (inherits(plot, "gglist")) {
    heights <- attr(plot, "params")$heights
    if (is.null(heights)) heights <- if (length(plot) == 3) c(1, 0.2, 0.6) else rep(1, length(plot))
    plot <- patchwork::wrap_plots(plot, ncol = 1, heights = heights)
  }
  if (inherits(plot, "patchwork")) return(patchwork::patchworkGrob(plot))
  if (inherits(plot, "ggplot")) return(ggplot2::ggplotGrob(plot))
  if (inherits(plot, "pheatmap")) return(plot$gtable)
  if (inherits(plot, c("grob", "gtable", "gTree"))) return(plot)
  stop("不支持的面板类型：", paste(class(plot), collapse = ", "))
}

apply_legacy_theme <- function(plot, style) {
  if (inherits(plot, "gglist")) {
    for (i in seq_along(plot)) plot[[i]] <- plot[[i]] + theme_project(style)
    return(plot)
  }
  if (inherits(plot, "patchwork")) return(plot & theme_project(style))
  plot + theme_project(style)
}

variable_palette <- function(style, variable, values) {
  palette <- style$variable_colors[[variable]]
  if (is.null(palette) && variable == "Condition") palette <- style$condition_colors
  if (is.null(palette)) return(NULL)
  palette <- unlist(palette)
  if (!all(unique(as.character(values)) %in% names(palette))) return(NULL)
  palette
}
