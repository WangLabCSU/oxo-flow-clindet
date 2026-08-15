#!/usr/bin/env Rscript
# smk.R — S4 shim that lets verbatim upstream clindet R scripts (which talk to
# Snakemake through the global `snakemake` S4 object) run under oxo-flow.
#
# Usage:
#   Rscript scripts/smk.R --script scripts/merge_maf.R \
#       --input maf1 file1 file2 ... --param dir DIR --output maf OUT
#
# Flags:
#   --script PATH   R script to run (sys.source'd verbatim)
#   --input NAME V...   named input slot (multiple values = vector)
#   --output NAME V...  named output slot
#   --param NAME V...   named params slot ("TRUE"/"FALSE" -> logical,
#                       integer strings -> integer)
#   --wildcard NAME V... named wildcard slot
#   --log NAME V...     named log slot
#   --threads N         threads slot
#
# Missing `--param` keys stay NULL, which the upstream scripts treat as
# "not provided" (e.g. the CNV/QC params of rmd.R).

args <- commandArgs(trailingOnly = TRUE)

is_flag <- function(x) startsWith(x, "--")

build_object <- function(args) {
  obj <- list(
    input = list(),
    output = list(),
    params = list(),
    wildcards = list(),
    log = list(),
    threads = 1L
  )
  script <- NULL
  i <- 1
  while (i <= length(args)) {
    flag <- sub("^--", "", args[i])
    if (flag == "script") {
      script <- args[i + 1]
      i <- i + 2
      next
    }
    if (flag == "threads") {
      obj$threads <- as.integer(args[i + 1])
      i <- i + 2
      next
    }
    # --input/--output/--param/--wildcard/--log NAME VALUES...
    if (flag %in% c("input", "output", "param", "wildcard", "log")) {
      name <- args[i + 1]
      vals <- character(0)
      j <- i + 2
      while (j <= length(args) && !is_flag(args[j])) {
        vals <- c(vals, args[j])
        j <- j + 1
      }
      converted <- lapply(vals, function(v) {
        if (v == "TRUE") return(TRUE)
        if (v == "FALSE") return(FALSE)
        if (grepl("^[0-9]+$", v)) return(as.integer(v))
        v
      })
      if (flag == "input") obj$input[[name]] <- if (length(converted) == 1) converted[[1]] else converted
      if (flag == "output") obj$output[[name]] <- if (length(converted) == 1) converted[[1]] else converted
      if (flag == "param") obj$params[[name]] <- if (length(converted) == 1) converted[[1]] else converted
      if (flag == "wildcard") obj$wildcards[[name]] <- if (length(converted) == 1) converted[[1]] else converted
      if (flag == "log") obj$log[[name]] <- if (length(converted) == 1) converted[[1]] else converted
      i <- j
      next
    }
    stop("smk.R: unknown flag --", flag)
  }
  if (is.null(script)) stop("smk.R: --script is required")
  list(obj = obj, script = script)
}

parsed <- build_object(args)

# Build the S4 class once; upstream scripts use `snakemake@input[['x']]`,
# `snakemake@params[['y']]`, `snakemake@output[['z']]`, `snakemake@log[[1]]`,
# `snakemake@wildcards.x`, `snakemake@threads`.
setClass(
  "snakemake_obj",
  representation(
    input = "list",
    output = "list",
    params = "list",
    wildcards = "list",
    log = "list",
    threads = "numeric"
  )
)
snakemake <- new(
  "snakemake_obj",
  input = parsed$obj$input,
  output = parsed$obj$output,
  params = parsed$obj$params,
  wildcards = parsed$obj$wildcards,
  log = parsed$obj$log,
  threads = parsed$obj$threads
)

# Run the verbatim upstream script in the global environment where `snakemake`
# is visible (upstream scripts reference it as a bare symbol).
sys.source(parsed$script, envir = globalenv())
