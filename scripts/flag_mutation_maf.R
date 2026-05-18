suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(GenomeInfoDb)
})

# -----------------------------
# Read arguments from Snakemake
# -----------------------------

maf_file <- snakemake@input[["maf"]]
bed_list_file <- snakemake@input[["bed_list"]]
output_file <- snakemake@output[["maf"]]

genome_version <- snakemake@params[["genome_version"]]
chr_col <- snakemake@params[["chr_col"]]
start_col <- snakemake@params[["start_col"]]
end_col <- snakemake@params[["end_col"]]
keep_chr_prefix <- snakemake@params[["keep_chr_prefix"]]
ignore_strand <- snakemake@params[["ignore_strand"]]
max_region_labels <- snakemake@params[["max_region_labels"]]

log_file <- snakemake@log[[1]]

# redirect messages to log
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

message("Input MAF: ", maf_file)
message("BED list: ", bed_list_file)
message("Output: ", output_file)
message("Genome version: ", genome_version)
message("chr_col: ", chr_col)
message("start_col: ", start_col)
message("end_col: ", end_col)
message("keep_chr_prefix: ", keep_chr_prefix)
message("ignore_strand: ", ignore_strand)

# -----------------------------
# Helper functions
# -----------------------------

safe_colname <- function(x) {
  x <- make.names(x)
  x <- gsub("\\.", "_", x)
  x
}

maf_to_granges <- function(maf,
                           chr_col = "Chromosome",
                           start_col = "Start_Position",
                           end_col = "End_Position",
                           keep_chr_prefix = TRUE) {
  required_cols <- c(chr_col, start_col, end_col)
  missing_cols <- setdiff(required_cols, colnames(maf))

  if (length(missing_cols) > 0) {
    stop("MAF missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  start <- as.integer(maf[[start_col]])
  end <- as.integer(maf[[end_col]])

  if (any(is.na(start))) {
    stop("Start position contains NA in column: ", start_col)
  }

  end[is.na(end) | end < start] <- start[is.na(end) | end < start]

  gr <- GRanges(
    seqnames = maf[[chr_col]],
    ranges = IRanges(start = start, end = end)
  )

  mcols(gr)$maf_index <- seq_len(nrow(maf))

  if (keep_chr_prefix) {
    seqlevelsStyle(gr) <- "UCSC"
  } else {
    seqlevelsStyle(gr) <- "NCBI"
  }

  gr
}

read_bed_as_granges <- function(bed_file,
                                annotation_name,
                                keep_chr_prefix = TRUE) {
  if (!file.exists(bed_file)) {
    stop("BED file does not exist: ", bed_file)
  }

  bed_gr <- rtracklayer::import(bed_file, format = "BED")
  if (keep_chr_prefix) {
    seqlevelsStyle(bed_gr) <- "UCSC"
  } else {
    seqlevelsStyle(bed_gr) <- "NCBI"
  }

  mcols(bed_gr)$bed_source <- annotation_name

  bed_gr
}

annotate_one_bed <- function(maf,
                             maf_gr,
                             bed_file,
                             annotation_name,
                             keep_chr_prefix = TRUE,
                             ignore_strand = TRUE,
                             max_region_labels = 5) {
  safe_name <- safe_colname(annotation_name)

  message("Annotating: ", annotation_name)
  message("BED: ", bed_file)

  bed_gr <- read_bed_as_granges(
    bed_file = bed_file,
    annotation_name = annotation_name,
    keep_chr_prefix = keep_chr_prefix
  )

  hits <- findOverlaps(
    query = maf_gr,
    subject = bed_gr,
    ignore.strand = ignore_strand
  )

  n <- length(maf_gr)

  flag <- rep(FALSE, n)
  overlap_count <- rep(0L, n)
  overlap_bases <- rep(0L, n)
  overlap_regions <- rep(NA_character_, n)

  if (length(hits) > 0) {
    qh <- queryHits(hits)
    sh <- subjectHits(hits)

    flag[unique(qh)] <- TRUE

    count_dt <- data.table(
      maf_index = qh
    )[
      ,
      .N,
      by = maf_index
    ]

    overlap_count[count_dt$maf_index] <- count_dt$N

    inter_gr <- pintersect(
      maf_gr[qh],
      bed_gr[sh],
      ignore.strand = ignore_strand
    )

    base_dt <- data.table(
      maf_index = qh,
      overlap_width = width(inter_gr)
    )[
      ,
      .(
        overlap_bases = sum(overlap_width, na.rm = TRUE)
      ),
      by = maf_index
    ]

    overlap_bases[base_dt$maf_index] <- base_dt$overlap_bases

    region_label <- paste0(
      as.character(seqnames(bed_gr[sh])),
      ":",
      start(bed_gr[sh]),
      "-",
      end(bed_gr[sh])
    )

    region_dt <- data.table(
      maf_index = qh,
      region = region_label
    )[
      ,
      .(
        regions = paste(
          unique(region)[seq_len(min(.N, max_region_labels))],
          collapse = ";"
        )
      ),
      by = maf_index
    ]

    overlap_regions[region_dt$maf_index] <- region_dt$regions
  }

  maf[[paste0("in_", safe_name)]] <- flag
  maf[[paste0(safe_name, "_overlap_count")]] <- overlap_count
  maf[[paste0(safe_name, "_overlap_bases")]] <- overlap_bases
  maf[[paste0(safe_name, "_overlap_regions")]] <- overlap_regions

  maf
}

run_region_annotation <- function(maf_file,
                                  bed_list_file,
                                  output_file,
                                  chr_col = "Chromosome",
                                  start_col = "Start_Position",
                                  end_col = "End_Position",
                                  keep_chr_prefix = TRUE,
                                  ignore_strand = TRUE,
                                  max_region_labels = 5) {
  maf <- fread(maf_file, data.table = TRUE)

  bed_list <- fread(bed_list_file, data.table = TRUE)

  required_bed_cols <- c("name", "path")
  missing_bed_cols <- setdiff(required_bed_cols, colnames(bed_list))

  if (length(missing_bed_cols) > 0) {
    stop("BED list must contain columns: name, path")
  }

  message("Number of variants in MAF: ", nrow(maf))
  message("Number of BED tracks: ", nrow(bed_list))

  maf_gr <- maf_to_granges(
    maf = maf,
    chr_col = chr_col,
    start_col = start_col,
    end_col = end_col,
    keep_chr_prefix = keep_chr_prefix
  )

  for (i in seq_len(nrow(bed_list))) {
    maf <- annotate_one_bed(
      maf = maf,
      maf_gr = maf_gr,
      bed_file = bed_list$path[i],
      annotation_name = bed_list$name[i],
      keep_chr_prefix = keep_chr_prefix,
      ignore_strand = ignore_strand,
      max_region_labels = max_region_labels
    )
  }

  in_cols <- paste0("in_", safe_colname(bed_list$name))
  in_cols <- intersect(in_cols, colnames(maf))

  maf[
    ,
    bed_annotation_count := rowSums(.SD, na.rm = TRUE),
    .SDcols = in_cols
  ]

  maf[
    ,
    bed_annotation_labels := apply(.SD, 1, function(x) {
      labs <- names(x)[as.logical(x)]
      labs <- sub("^in_", "", labs)

      if (length(labs) == 0) {
        return("NONE")
      }

      paste(labs, collapse = ";")
    }),
    .SDcols = in_cols
  ]

  outdir <- dirname(output_file)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  fwrite(
    maf,
    output_file,
    sep = "\t",
    quote = FALSE,
    na = ""
  )

  message("Done.")
  message("Output written to: ", output_file)

  invisible(maf)
}

# -----------------------------
# Run
# -----------------------------

run_region_annotation(
  maf_file = maf_file,
  bed_list_file = bed_list_file,
  output_file = output_file,
  chr_col = chr_col,
  start_col = start_col,
  end_col = end_col,
  keep_chr_prefix = keep_chr_prefix,
  ignore_strand = ignore_strand,
  max_region_labels = max_region_labels
)