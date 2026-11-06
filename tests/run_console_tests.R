Sys.setenv(CONSOLE_PROJECT_ROOT = normalizePath(getwd()))
testthat::test_dir("tests/console", reporter = "summary", stop_on_failure = TRUE)
