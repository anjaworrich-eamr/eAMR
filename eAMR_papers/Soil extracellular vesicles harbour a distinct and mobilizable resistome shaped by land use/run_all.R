# =====================================================================
#  run_all.R - run all scripts in order
# ---------------------------------------------------------------------
#  Each script runs in its own R process; a failing script does not stop the others.
#  Logs: logs/<script>.log; overview: logs/run_all_status.tsv
#  Use: open soil-EV-resistome.Rproj, then source("scripts/run_all.R")
# =====================================================================
if (!file.exists("soil-EV-resistome.Rproj")) stop("Open soil-EV-resistome.Rproj before running this script")
dir.create("logs", showWarnings = FALSE)
scripts <- c("00_prepare_inputs",
             "Fig1", "Fig2", "Fig3", "Fig4", "Fig5",
             "FigS2", "FigS3", "FigS4", "FigS5", "FigS6", "FigS7", "FigS8", "FigS9", "FigS10_S11", "FigS12_S13",
             "TableS1", "TableS2", "TableS3")
rscript <- file.path(R.home("bin"), "Rscript")
status  <- file.path("logs", "run_all_status.tsv")
cat("script\tstatus\tseconds\n", file = status)
for (s in scripts) {
  t0  <- Sys.time()
  log <- file.path("logs", paste0(s, ".log"))
  code <- system2(rscript, c("-e", shQuote(sprintf("options(repos = c(CRAN = 'https://cloud.r-project.org')); source('scripts/%s.R')", s))),
                  stdout = log, stderr = log)
  cat(sprintf("%s\t%s\t%.0f\n", s, ifelse(code == 0, "ok", "FAILED"), as.numeric(difftime(Sys.time(), t0, units = "secs"))),
      file = status, append = TRUE)
}
