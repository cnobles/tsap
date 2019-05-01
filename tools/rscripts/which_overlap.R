#' Script to generate summary report for VIVI analysis
#' Required packages:

# Initial options and packages -------------------------------------------------
options(stringsAsFactors = FALSE, scipen = 99)
packs <- c("argparse", "magrittr", "pander")
loaded_packs <- suppressMessages(sapply(packs, require, character.only = TRUE))
if(any(!loaded_packs)){
  cat("Check below for required packages:")
  print(loaded_packs)
  stop("Check dependancies.")
}
panderOptions("table.style", "simple")
panderOptions("table.split.table", Inf)


# Argument parser --------------------------------------------------------------
parser <- ArgumentParser(
  description = "R-based tool to identify overlapping reads based on reference alignment.")
parser$add_argument(
  "-f", "--fwd.algns", nargs = 1, type = "character", 
  help = "Alignments from the forward read to references. BAM file.")
parser$add_argument(
  "-r", "--rev.algns", nargs = 1, type = "character", 
  help = "Alignments from the reverse read to reference.")
parser$add_argument(
  "--fi", nargs = 1, type = "character", 
  help = "Output path for report. No extension required.")
parser$add_argument(
  "--ri", nargs = 1, type = "character",
  help = "Run specific config file(s) in yaml format.")
parser$add_argument(
  "-o", "--output", nargs = 1, type = "character",
  help = "Output file with readnames that overlap between foward and reverse reads.")
parser$add_argument(
  "--ref", nargs = 1, type = "character",
  help = "Reference sequences specific to panel.")
parser$add_argument(
  "-n", "--min", nargs = 1, type = "integer", default = 10,
  help = "Minimum overlap between foward and reverse alignments. Default 10.")

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Print out inputs for user / logs
input_table <- data.frame(
  "Variable" = names(args), 
  "Value" = sapply(
    seq_along(args), function(i) paste(args[[i]], collapse = ", ")))
input_table <- input_table[
  match(
    c("fwd.algns", "fi", "rev.algns", "ri", 
      "output", "ref", "min"), 
    input_table$Variable),]
input_table$Variable <- paste0(format(input_table$Variable), " :")
pander("Reference based overlap:\n")
pander(
  data.frame(input_table, row.names = NULL),
  row.names = FALSE, justify = "left")

# Load additional supporting packages ------------------------------------------
message("Loading dependencies.")
add_packs <- c(
  "argparse", "magrittr", "pander", "data.table", "Biostrings", "GenomicRanges",
  "Rsamtools", "tidyverse")
add_packs <- c(packs, add_packs)

add_packs_loaded <- suppressMessages(
  sapply(add_packs, require, character.only = TRUE))
if(!all(add_packs_loaded)){
  pandoc.table(data.frame(
    "R-Packages" = names(add_packs_loaded), 
    "Loaded" = add_packs_loaded, 
    row.names = NULL))
  stop("Check dependancies.")
}

# Supporting functions ---------------------------------------------------------
bam_params <- c("qname", "flag", "rname", "pos", "qwidth", "mapq", "cigar")

cigarRanges <- function(rname, cigar, pos, sym, seq.info){
  sym <- unlist(strsplit(sym, ""))
  
  # Parse cigar strings
  comp_sym <- IRanges::CharacterList(stringr::str_split(cigar, "[0-9]+"))
  comp_sym <- comp_sym[which(comp_sym != "")]
  comp_val <- IRanges::IntegerList(stringr::str_split(cigar, "[A-z]"))
  comp_val <- comp_val[which(!is.na(comp_val))]
  
  # Remove soft and hard clipping
  keep_log <- !comp_sym %in% c("H", "S")
  comp_sym <- comp_sym[keep_log]
  comp_val <- comp_val[keep_log]
  
  # Extract insertions
  insert_log <- comp_sym == "I"
  insert_val <- comp_val[insert_log]
  insert_sym <- comp_sym[insert_log]
  mdi_val <- IRanges::ifelse2(!insert_log, comp_val, 0)
  
  # Calculate starts and ends
  sum_val <- cumsum(mdi_val)
  mdi_start <- sum_val
  mdi_start@unlistData <- c(
    0, mdi_start@unlistData[seq_len(length(mdi_start@unlistData)-1)] + 1)
  mdi_start@unlistData[start(mdi_start@partitioning)] <- pos
  mdi_end <- sum_val + pos - 1
  mdi_end <- IRanges::ifelse2(insert_log, mdi_end + 1, mdi_end)
  
  # Construct GRanges for M and D
  gr_md <- GenomicRanges::GRanges(
    seqnames = rep(rname, lengths(mdi_start[!insert_log])),
    ranges = IRanges::IRanges(
      start = unlist(mdi_start[!insert_log]), 
      end = unlist(mdi_end[!insert_log])),
    strand = "+",
    seqinfo = seq.info,
    "index" = rep(seq_along(rname), lengths(mdi_start[!insert_log])),
    "symbol" = unlist(comp_sym[!insert_log]),
    "in.len" = NA)
  
  # Construct GRanges for I
  gr_i <- GenomicRanges::GRanges(
    seqnames = rep(rname, lengths(mdi_start[insert_log])),
    ranges = IRanges::IRanges(
      start = unlist(mdi_start[insert_log]), 
      end = unlist(mdi_end[insert_log])),
    strand = rep("+", sum(lengths(mdi_start[insert_log]))),
    seqinfo = seq.info,
    "index" = rep(seq_along(rname), lengths(mdi_start[insert_log])),
    "symbol" = unlist(comp_sym[insert_log]),
    "in.len" = unlist(insert_val))
  
  # Combine GRanges
  gr <- c(gr_md, gr_i)
  
  # Filter output
  gr <- gr[gr$symbol %in% sym]
  
  # Return GRangesList
  split(gr, gr$index)
}

import_bam <- function(bam.path, bai.path, params){
  bam_env <- new.env()
  bam_env$bam <- unlist(
    Rsamtools::scanBam(
      bam.path, index = bai.path, 
      param = Rsamtools::ScanBamParam(what = params)),
    recursive = FALSE)
  bam_env$bam_df <- as.data.frame(bam_env$bam[seq_along(params)])
  bam <- bam_env$bam_df
  names(bam) <- tolower(names(bam))
  null_algns <- bam[which(is.na(bam$rname)),]
  bam[!is.na(bam$rname),]
}

# Load references --------------------------------------------------------------
ref_seqs <- Biostrings::readDNAStringSet(args$ref, format = "fasta")
ref_seqinfo <- GenomeInfoDb::Seqinfo(
  seqnames = names(ref_seqs), 
  seqlengths = Biostrings::width(ref_seqs), 
  isCircular = rep(FALSE, length(ref_seqs)))

# Import alignment files -------------------------------------------------------
fwd.bam <- import_bam(args$fwd.algns, args$fi, bam_params)
fwd.cig <- with(
  fwd.bam, cigarRanges(rname, cigar, pos, sym = "MD", seq.info = ref_seqinfo))
names(fwd.cig) <- fwd.bam$qname
rev.bam <- import_bam(args$rev.algns, args$ri, bam_params)
rev.cig <- with(
  rev.bam, cigarRanges(rname, cigar, pos, sym = "MD", seq.info = ref_seqinfo))
names(rev.cig) <- rev.bam$qname

# Intersecting name and ranges -------------------------------------------------
common_names <- intersect(names(fwd.cig), names(rev.cig))
fwd.cig <- GenomicRanges::reduce(fwd.cig[common_names])
rev.cig <- GenomicRanges::reduce(rev.cig[common_names])

idx <- which(max(end(fwd.cig)) - min(start(rev.cig)) >= args$min)

ovlp_reads <- names(idx)

# Write names to text file -----------------------------------------------------
write.csv(
  data.frame(ovlp_reads), 
  file = args$output, 
  row.names = FALSE, quote = FALSE)

q()

