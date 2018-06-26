#' Script for processing alignments
#' Requires: 
#'   argparse, Rsamtools, magrittr, tidyverse, data.table, yaml, pander, 
#'   GenomeInfoDb, Biostrings, BSgenome, GenomicRanges, IRanges

packs <- c("argparse", "magrittr", "pander")
loaded_packs <- suppressMessages(sapply(packs, require, character.only = TRUE))
if(any(!loaded_packs)){
  cat("Check below for required packages:")
  print(loaded_packs)
  stop()
}
options(stringsAsFactors = FALSE, scipen = 99)
panderOptions("table.split.table", Inf)

# Argument parser --------------------------------------------------------------
parser <- ArgumentParser(
  description = "R-based post alignment analysis focused on identifying indels.")
parser$add_argument(
  "bamFile", type = "character", 
  help = "Path to bam file input.")
parser$add_argument(
  "-i", "--index", type = "character", 
  help = "Path to index of bam file input.")
parser$add_argument(
  "-k", "--key", type = "character", 
  help = "Path to key file for appending read count data.")
parser$add_argument(
  "-o", "--output", type = "character", 
  help = "Output path for unique alignments, csv format.")
parser$add_argument(
  "-a", "--chimera", type = "character", 
  help = "Output path for chimeric or translocation alignments.")
parser$add_argument(
  "-c", "--config", type = "character", 
  help = "Config file for run.")
parser$add_argument(
  "-r", "--ref", type = "character", 
  help = "Reference sequences specific to panel.")
parser$add_argument(
  "-t", "--target", type = "character", 
  help = "Panel target file, csv format.")

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

input_table <- data.frame(
  "Variable" = names(args), 
  "Value" = sapply(
    seq_along(args), function(i) paste(args[[i]], collapse = ", ")))
input_table <- input_table[
  match(
    c("bamFile", "index", "key", "output", "chimera", "config", "ref"), 
    input_table$Variable),]
input_table$Variable <- paste0(format(input_table$Variable), " :")
pander("Post-alignment Processing Inputs\n")
pander(
  data.frame(input_table, row.names = NULL),
  row.names = FALSE, justify = "left")


# Script parameters and additional functions -----------------------------------
bam_params <- c("qname", "flag", "rname", "pos", "qwidth", "mapq", "cigar")
bam_tags <- c("MD")

getGenome <- function(ref){
  gen_name <- grep(ref, unique(BSgenome::installed.genomes()), value = TRUE)
  
  if(length(gen_name) < 1){
    stop("No matched BSgenome installed. Please install.")
  }else if(length(gen_name) > 1){
    message("Installed matching genomes:\n")
    message(paste(gen_name, collapse = ", "))
    stop("Ambiguous match to requested genome. Please specify.")
  }
  
  suppressMessages(library(gen_name, character.only = TRUE))
  BiocGenerics::get(gen_name)
}

cntSym <- function(cigar, sym, max.only = FALSE){
  match_sym <- paste0("[0-9]+", sym)
  ex_mat <- stringr::str_extract_all(cigar, match_sym, simplify = TRUE)
  ex_mat <- as.numeric(gsub(sym, "", ifelse(nchar(ex_mat) == 0, 0, ex_mat)))
  ex_mat <- matrix(ex_mat, nrow = length(cigar))
  if(ncol(ex_mat) == 0) ex_mat <- matrix(rep(0, length(cigar)), ncol = 1)
  
  if(max.only){
    return(ex_mat[
      matrix(c(seq_len(nrow(ex_mat)), max.col(ex_mat, "first")), ncol = 2)])
  }else{
    return(rowSums(ex_mat, na.rm = TRUE))
  }
}

getFlankingSeqs <- function(posid, ref, flk = 30){
  seq <- stringr::str_extract(posid, "[\\w]+")
  pos <- as.numeric(stringr::str_extract(posid, "[0-9]+$"))
  gr <- GenomicRanges::GRanges(
    seqnames = seq, 
    ranges = IRanges::IRanges(start = pos, end = pos + flk - 1),
    strand = "+",
    seqinfo = BSgenome::seqinfo(ref))
  BSgenome::getSeq(ref, gr)
}

mindexToGranges <- function(mindex, strand, ref = NULL){
  if(is.null(mindex@NAMES)) stop("NAMES column not found for seqnames.")
  ir <- unlist(mindex)
  if(is.null(ref)){
    return(unname(GenomicRanges::GRanges(
      seqnames = names(ir),
      ranges = ir,
      strand = rep(strand, length(ir)))))
  }else{
    return(unname(GenomicRanges::GRanges(
      seqnames = names(ir),
      ranges = ir,
      strand = rep(strand, length(ir)),
      seqinfo = ref)))
  }
}

normalizeTargets <- function(posid, base.ref, tar.seqs, tar.seqinfo,
                              flk = 30, return = "posid"){
  stopifnot(return %in% c("posid", "granges"))
  tar_oligos <- getFlankingSeqs(posid, base.ref, flk = flk)
  names(tar_oligos) <- names(posid)
  tar_matches <- unlist(GenomicRanges::GenomicRangesList(lapply(
    seq_along(tar_oligos), function(i){
      x <- Biostrings::vmatchPattern(
        pattern = as.character(tar_oligos[i]), 
        subject = tar.seqs[match(names(tar_oligos)[i], names(tar.seqs))])
      mindexToGranges(x, "+", tar.seqinfo)
    })))
  tar_pos <- flank(tar_matches, width = -1, start = TRUE)
  if(return == "posid"){
    return(structure(
      paste0(
        GenomicRanges::seqnames(tar_pos), ":", GenomicRanges::start(tar_pos)),
      names = names(posid)))
  }else if(return == "granges"){
    return(tar_pos)
  }
}

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
  mdi_val <- ifelse(!insert_log, comp_val, 0)
  
  # Calculate starts and ends
  sum_val <- cumsum(mdi_val)
  mdi_start <- sum_val
  mdi_start@unlistData <- c(
    0, mdi_start@unlistData[seq_len(length(mdi_start@unlistData)-1)] + 1)
  mdi_start@unlistData[start(mdi_start@partitioning)] <- pos
  mdi_end <- sum_val + pos - 1
  mdi_end <- ifelse(insert_log, mdi_end + 1, mdi_end)
  
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

# Import inputs for analysis ---------------------------------------------------
## Load config file
config <- yaml::yaml.load_file(args$config)

## Sequence Reference
ref_seqs <- Biostrings::readDNAStringSet(args$ref, format = "fasta")
ref_seqinfo <- GenomeInfoDb::Seqinfo(
  seqnames = names(ref_seqs), 
  seqlengths = Biostrings::width(ref_seqs), 
  isCircular = rep(FALSE, length(ref_seqs)), 
  genome = rep(config$RefGenome, length(ref_seqs)))

base_gen <- getGenome(config$RefGenome)

## Panel targets
panel_targets <- read.csv(
  file.path(config$Install_Directory, config$Panel_Path))
targets <- structure(panel_targets$locus, names = panel_targets$target)
norm_targets <- normalizeTargets(
  targets, base_gen, ref_seqs, ref_seqinfo, return = "granges")

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
null_algns <- bam[which(is.na(bam$rname)),]
bam <- bam[!is.na(bam$rname),]

# Summarize data and conduct prelimininary analysis ----------------------------
## Generate data.frame of unique alignment data
uniq_algns <- dplyr::group_by(bam, qname) %>%
  dplyr::filter(n() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    in.cnt = cntSym(cigar, "I"),
    in.max = cntSym(cigar, "I", TRUE),
    del.cnt = cntSym(cigar, "D"),
    del.max = cntSym(cigar, "D", TRUE))

uniq_indel_grl <- uniq_algns %$%
  cigarRanges(
    rname, cigar, pos, sym = c("D", "I"), seq.info = ref_seqinfo)

indel_algns <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(
    norm_targets, uniq_indel_grl, maxgap = config$maxDistFromEdit))

uniq_algns <- dplyr::mutate(
  uniq_algns, tar.indel = seq_len(n()) %in% 
    as.numeric(names(uniq_indel_grl[indel_algns])))

if(!is.null(args$key)){
  uniq_algns <- dplyr::left_join(uniq_algns, key, by = c("qname" = "seqID"))
}

## Generate data.frame of chimeric or translocation alignments
if(!is.null(args$chimera)){
  chim_algns <- dplyr::group_by(bam, qname) %>%
    dplyr::filter(n() > 1, n_distinct(rname) > 1) %>%
    dplyr::arrange(pos) %>%
    dplyr::summarise(
      pri.seq = dplyr::first(rname[which(pos == min(pos))]),
      pri.start = dplyr::first(pos[which(pos == min(pos))]),
      pri.width = dplyr::first(cntSym(cigar, "M")[which(pos == min(pos))]),
      sec.seq = dplyr::last(rname[which(pos == max(pos))]),
      sec.start = dplyr::last(pos[which(pos == max(pos))]),
      sec.width = dplyr::last(cntSym(cigar, "M")[which(pos == max(pos))])) %>%
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