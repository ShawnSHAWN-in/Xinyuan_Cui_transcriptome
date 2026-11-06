Sys.setenv(CONSOLE_PROJECT_ROOT = normalizePath(getwd()))
loadNamespace("targets")
for (package in c("crew", "visNetwork", "parallelly")) loadNamespace(package)
testthat::test_dir("tests/console", reporter = "summary", stop_on_failure = TRUE)
