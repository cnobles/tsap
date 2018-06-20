# Script for processing alignments
# Requires: argparse, Rsamtools, magrittr, tidyverse, data.table
options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))
suppressMessages(library("magrittr"))

# Argument parser --------------------------------------------------------------
parser <- ArgumentParser(
  description = "R-based post alignment analysis focused on identifying indels.")
parser$add_argument(
  "bamFile", type = "character", help = "Path to bam file input.")
parser$add_argument(
  "-i", "--index", type = "character", help = "Path to index of bam file input.")
parser$add_argument(
  "-k", "--key", type = "character", help = "Path to key file for appending read count data.")
parser$add_argument(
  "-o", "--output", type = "character", help = "Output path for unique alignments, csv format.")
parser$add_argument(
  "-c", "--chimera", type = "character", help = "Output path for chimeric or translocation alignments.")

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

input_table <- data.frame(
  "Variable" = names(args), 
  "Value" = sapply(
    seq_along(args), function(i) paste(args[[i]], collapse = ", ")))
input_table <- input_table[
  match(
    c("bamFile", "index", "key", "output", "chimera"), 
    input_table$Variable),]
input_table$Variable <- paste0(format(input_table$Variable), " :")
cat("Post-alignment Processing Inputs\n")
print.data.frame(
  data.frame(input_table, row.names = NULL),
  row.names = FALSE, right = FALSE)


# Script parameters and additional functions -----------------------------------
bam_params <- c("qname", "flag", "rname", "pos", "qwidth", "mapq", "cigar")
bam_tags <- c("MD")

cigar_match <- function(cigar){
  rowSums(matrix(as.numeric(gsub("M", "", stringr::str_extract_all(
    cigar, "[0-9]+M", simplify = TRUE))), 
    nrow = length(cigar)), na.rm = TRUE)
}

insert_cnt <- function(cigar){
  rowSums(matrix(as.numeric(gsub("I", "", stringr::str_extract_all(
    cigar, "[0-9]+I", simplify = TRUE))), 
    nrow = length(cigar)), na.rm = TRUE)
}

delete_cnt <- function(cigar){
  rowSums(matrix(as.numeric(gsub("D", "", stringr::str_extract_all(
    cigar, "[0-9]+D", simplify = TRUE))), 
    nrow = length(cigar)), na.rm = TRUE)
}

# Import inputs for analysis ---------------------------------------------------
## Load key file
if(!is.null(args$key)){
  key <- data.table::fread(args$key, header = TRUE, data.table = FALSE) %>%
    dplyr::group_by(seqID) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::ungroup()
}

## Load bam file and collect relevent data
bam_env <- new.env()
bam_env$bam <- unlist(
  Rsamtools::scanBam(
    args$bamFile, 
    index = args$index, 
    param = Rsamtools::ScanBamParam(what = bam_params, tag = bam_tags)),
  recursive = FALSE) 
bam_env$bam_df <- as.data.frame(bam_env$bam[seq_along(bam_params)])
null <- lapply(bam_tags, function(t, src){
  bam_env$bam_df[,t] <- src$tag[[t]]
  }, src = bam_env$bam)
bam <- bam_env$bam_df
names(bam) <- tolower(names(bam))

## Generate data.frame of unique alignment data
uniq_algns <- dplyr::group_by(bam, qname) %>%
  dplyr::filter(n() == 1) %>%
  dplyr::mutate(
    in.cnt = insert_cnt(cigar),
    del.cnt = delete_cnt(cigar))

if(!is.null(args$key)){
  uniq_algns <- dplyr::left_join(uniq_algns, key, by = c("qname" = "seqID"))
}

## Generate data.frame of chimeric or translocation alignments
if(!is.null(args$chimera)){
  chim_algns <- dplyr::group_by(bam, qname) %>%
    dplyr::filter(n() > 1, dplyr::n_distinct(rname) > 1) %>%
    dplyr::arrange(pos) %>%
    dplyr::summarise(
      pri.seq = rname[which(pos == min(pos))],
      pri.start = pos[which(pos == min(pos))],
      pri.width = cigar_match(cigar)[which(pos == min(pos))],
      sec.seq = rname[which(pos == max(pos))],
      sec.start = pos[which(pos == max(pos))],
      sec.width = cigar_match(cigar)[which(pos == max(pos))]) %>%
    dplyr::ungroup()
  if(!is.null(args$key)){
    chim_algns <- dplyr::left_join(chim_algns, key, by = c("qname" = "seqID"))
  }
}

# Write output files -----------------------------------------------------------
write.csv(uniq_algns, file = args$output, quote = FALSE, row.names = FALSE)

if(!is.null(args$chimera)){
  write.csv(chim_algns, file = args$chimera, quote = FALSE, row.names = FALSE)
}

q()