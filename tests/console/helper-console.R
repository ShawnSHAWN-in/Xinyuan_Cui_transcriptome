project_root <- normalizePath(Sys.getenv("CONSOLE_PROJECT_ROOT", unset = "../.."))
source(file.path(project_root, "console", "backend.R"))
source(file.path(project_root, "console", "app.R"))
source(file.path(project_root, "scripts", "console.R"))

console_fixture <- function() {
  root <- tempfile("console_中文 space_'_")
  dir.create(root)
  file.copy(file.path(project_root, c("config", "R", "console", "scripts")), root, recursive = TRUE)
  grid <- expand.grid(Replicate = 1:3, Condition = c("Cat", "Leopard", "Tiger", "NC"), Time_f = c("06", "12", "24"), stringsAsFactors = FALSE)
  sample <- paste(grid$Condition, grid$Time_f, grid$Replicate, sep = "_")
  metadata <- data.frame(Sample = sample, Group = paste(grid$Condition, grid$Time_f, sep = "_"))
  data.table::fwrite(metadata, file.path(root, "samples.tsv"), sep = "\t")
  counts <- data.frame(Geneid = "gene1", Length = 1000, check.names = FALSE)
  for (name in sample) counts[[name]] <- 10
  data.table::fwrite(counts, file.path(root, "counts.tsv"), sep = "\t")
  writeLines("annotation fixture", file.path(root, "reference.gtf"))
  writeLines("annotation fixture", file.path(root, "gene2go.tsv"))
  paths <- list(count_file = "counts.tsv", sample_file = "samples.tsv", gtf_file = "reference.gtf", gene2go_file = "gene2go.tsv", output_dir = "output")
  yaml::write_yaml(paths, file.path(root, "config", "paths.yml"))
  writeLines(c("list(targets::tar_target(qc_files, {",
    "paths <- yaml::read_yaml(file.path(Sys.getenv('TRANSCRIPTOME_CONFIG_DIR'), 'paths.yml'))",
    "out <- file.path(Sys.getenv('TRANSCRIPTOME_PROJECT_DIR'), paths$output_dir)",
    "dir.create(out, recursive = TRUE, showWarnings = FALSE)",
    "writeLines('fixture result', file.path(out, 'done.tsv'))",
    "file.path(out, 'done.tsv')", "}, format = 'file'))"), file.path(root, "_targets.R"))
  context <- console_context_lst(root)
  list(root = root, context = context, config = console_effective_config_lst(context))
}

wait_console <- function(controller, timeout = 30) {
  deadline <- Sys.time() + timeout
  repeat {
    state <- console_poll_lst(controller)
    if (state$state != "RUNNING") return(state)
    if (Sys.time() > deadline) {
      console_stop_run(controller)
      stop("后台验证超时")
    }
    Sys.sleep(0.05)
  }
}
