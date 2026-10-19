test_root <- Sys.getenv("ENGINEERING_PROJECT_ROOT")
fixture_path <- function(name) file.path(test_root, "tests", "fixtures", name)
config <- function(name) read_config_lst(file.path(test_root, "config", paste0(name, ".yml")))
legacy_files <- file.path(test_root, "analysis_scripts", "lib", c("Rscript_R6_transcriptome_base.R", "Rscript_R6_transcriptome_class.R",
                                                                "compatibility.R", "project_helpers.R"))
