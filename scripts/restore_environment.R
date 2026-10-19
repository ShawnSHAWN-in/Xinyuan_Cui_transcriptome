if (!requireNamespace("renv", quietly = TRUE)) stop("请先install.packages('renv')")
renv::restore(project = getwd(), prompt = FALSE)
message("环境已恢复；新R会话先运行renv::load()，再启动targets")
