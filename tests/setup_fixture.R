setup_targets_fixture <- function(project, root, quick = FALSE) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(project, c("R", "config", "analysis_scripts")), root, recursive = TRUE)
  file.copy(file.path(project, "_targets.R"), root)
  fixtures <- file.path(project, "tests", "fixtures")
  paths <- list(count_file = file.path(fixtures, "counts.tsv"), sample_file = file.path(fixtures, "samples.tsv"),
                gtf_file = file.path(fixtures, "reference.gtf"), gene2go_file = file.path(fixtures, "gene2go.tsv"), output_dir = "output")
  yaml::write_yaml(paths, file.path(root, "config", "paths.yml"))
  reference <- yaml::read_yaml(file.path(root, "config", "reference.yml"))
  reference$go_ancestors <- FALSE
  reference$kegg_term2gene <- file.path(fixtures, "KEGG_TERM2GENE.tsv")
  reference$kegg_term2name <- file.path(fixtures, "KEGG_TERM2NAME.tsv")
  yaml::write_yaml(reference, file.path(root, "config", "reference.yml"))
  figures <- yaml::read_yaml(file.path(root, "config", "figures.yml"))
  figures$export$dpi <- 36
  figures$export$max_pixels <- 600000
  figures$display$curve_top_each <- 1L
  figures$display$top_terms <- 3L
  figures$display$compare_max_terms <- 6L
  figures$display$top_genes_each <- 2L
  figures$display$enrichment_networks <- !quick
  figures$display$enrichment_upset <- !quick
  yaml::write_yaml(figures, file.path(root, "config", "figures.yml"))
  normalizePath(root)
}
