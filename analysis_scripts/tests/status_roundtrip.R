# 验证日志重复写入
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- normalizePath(file.path(dirname(sub("^--file=", "", arg[1])), ".."))
cfg <- list(output_dir = tempfile("status_roundtrip_"))
source(file.path(script_dir, "lib", "project_helpers.R"))
detail <- 'Package "ggupset" is required; quoted "text" stays unchanged'
for (i in seq_len(100)) record_status("test", paste0("item_", i), "FAILED", detail)
path <- file.path(cfg$output_dir, "10_objects_and_summary", "module_status.tsv")
tbl <- data.table::fread(path, quote = "", colClasses = "character", data.table = FALSE)
stopifnot(nrow(tbl) == 100L, all(tbl$Detail == detail), file.info(path)$size < 50000)
record_status("test", "item_1", "SUCCESS", "Recovered")
tbl <- data.table::fread(path, quote = "", colClasses = "character", data.table = FALSE)
stopifnot(nrow(tbl) == 100L, tbl$Status[tbl$Item == "item_1"] == "SUCCESS")
cat("STATUS_ROUNDTRIP_PASS\n")
