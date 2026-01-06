# Set R temp directory to repository .tmp/R instead of system /tmp
if (Sys.getenv("TMPDIR") == "") {
  repo_root <- normalizePath(".")
  tmp_dir <- file.path(repo_root, ".tmp", "R")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(TMPDIR = tmp_dir)
  Sys.setenv(TEMP = tmp_dir)
  Sys.setenv(TMP = tmp_dir)
}
