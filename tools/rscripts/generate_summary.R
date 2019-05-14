#' Script to generate summary report for TsAP analysis
#' Required packages:

# Initial options and packages -------------------------------------------------
options(stringsAsFactors = FALSE, scipen = 99)
packs <- c("argparse", "magrittr", "pander")

loaded_packs <- suppressMessages(sapply(packs, require, character.only = TRUE))

if( any(!loaded_packs) ){
  cat("Check below for required packages:")
  print(loaded_packs)
  stop("Check dependancies.")
}

panderOptions("table.style", "simple")
panderOptions("table.split.table", Inf)


# Argument parser --------------------------------------------------------------
parser <- ArgumentParser(
  description = "R-based post alignment analysis focused on identifying indels.")

parser$add_argument(
  "-u", "--unique", type = "character", 
  help = "Unique alignments output from processing.")

parser$add_argument(
  "-a", "--chimera", type = "character", 
  help = "Chimeric alignments output from processing.")

parser$add_argument(
  "-o", "--output", type = "character", 
  help = "Output path for report. No extension required.")

parser$add_argument(
  "-c", "--config", nargs = 1, type = "character",
  help = "Run specific config file(s) in yaml format.")

parser$add_argument(
  "-r", "--ref", nargs = 1, type = "character",
  help = "Reference sequences specific to panel.")

parser$add_argument(
  "-s", "--support", nargs = 1, type = "character",
  help = "Supplementary data input, csv or tsv format.")

parser$add_argument(
  "-f", "--figures", action = "store_true",
  help = "Generate figures along with output report (pdf and png formats).")

parser$add_argument(
  "-b", "--tables", action = "store_true",
  help = "Generate tables along with output report (csv formats)."
)

parser$add_argument(
  "-d", "--data", action = "store_true",
  help = "Data to generate the report will be saved as an R image with output.")

parser$add_argument(
  "-t", "--format", nargs = 1, type = "character", default = "html",
  help = "Output format for report. Either 'pdf' or 'html' (default).")

parser$add_argument(
  "--template", nargs = 1, type = "character", 
  default = "tools/rscripts/report_templates/report_template.Rmd",
  help = "File path to standard or custom iGUIDE report template."
)

parser$add_argument(
  "--install_dir", nargs = 1, type = "character", default = "TSAP_DIR",
  help = "VivI install directory path, do not change for normal applications."
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if( !dir.exists(args$install_dir) ){
  root_dir <- Sys.getenv(args$install_dir)
}else{
  root_dir <- args$install_dir
}

if( !dir.exists(root_dir) ){
  stop(paste0("\n  Cannot find install path to VivI: ", root_dir, ".\n"))
}else{
  args$install_dir <- root_dir
}

# Check for report output format
report_formats <- c("html" = "html_document", "pdf" = "pdf_document")

if( !args$format %in% names(report_formats) ){
  stop(
    "\n  Please input either 'html' or 'pdf' for format.\n",
    "\n  Other formats not currenlty supported.\n")
}

output_format <- report_formats[args$format]

# Resolve template file path.
if( file.exists(file.path(root_dir, args$template)) ){
  
  template_path <- normalizePath(file.path(root_dir, args$template))
  
}else if( file.exists(file.path(args$template)) ){
  
  template_path <- normalizePath(file.path(args$template))
  
}else{
  
  stop("\n  Cannot find template file: ", args$template, ".\n")
  
}

# Print out inputs for user / logs
input_table <- data.frame(
  "Variable" = names(args), 
  "Value" = sapply(
    seq_along(args), function(i) paste(args[[i]], collapse = ", "))
)

input_table <- input_table[
  match(
    c(
      "unique", "chimera", "output", "format", "config", "ref", "support", 
      "figures", "tables", "data", "template", "install_dir"
    ), 
    input_table$Variable),
]

input_table$Variable <- paste0(format(input_table$Variable), " :")

pander("Generate Summary Report Inputs:\n")

pander(
  data.frame(input_table, row.names = NULL),
  row.names = FALSE, justify = "left"
)


# Load remaining dependencies --------------------------------------------------
message("Loading dependencies.")

add_packs <- c(
  "magrittr", "tidyverse", "data.table", "Biostrings", "GenomicRanges", "knitr"
)

add_packs <- c(packs, add_packs)

add_packs_loaded <- suppressMessages(
  sapply(add_packs, require, character.only = TRUE)
)

if( !all(add_packs_loaded) ){

  pandoc.table(data.frame(
    "R-Packages" = names(add_packs_loaded), 
    "Loaded" = add_packs_loaded, 
    row.names = NULL))
  
  stop("Check dependancies.")
  
}

code_dir <- dirname(sub(
  "--file=", "", 
  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
)


# Supporting functions ---------------------------------------------------------
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
  
  ex_mat <- as.numeric(
    gsub(sym, "", IRanges::ifelse2(nchar(ex_mat) == 0, 0, ex_mat))
  )
  
  ex_mat <- matrix(ex_mat, nrow = length(cigar))
  
  if( ncol(ex_mat) == 0 ) ex_mat <- matrix(rep(0, length(cigar)), ncol = 1)
  
  if( max.only ){
    
    return(ex_mat[
      matrix(c(seq_len(nrow(ex_mat)), max.col(ex_mat, "first")), ncol = 2)
    ])
    
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
  mdi_val <- IRanges::ifelse2(!insert_log, comp_val, 0)
  
  # Calculate starts and ends
  sum_val <- cumsum(mdi_val) + pos - 1
  mdi_start <- sum_val
  mdi_start@unlistData <- c(
    0, mdi_start@unlistData[seq_len(length(mdi_start@unlistData)-1)] + 1)
  mdi_start@unlistData[start(mdi_start@partitioning)] <- pos
  mdi_end <- sum_val
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
  return(split(gr, gr$index))
}

pNums <- function(x, round = 3, digits = 3, big.mark = ",", ...){
  x <- ifelse(is.na(x), 0, x)
  if(x >= 100){
    return(format(x, big.mark = big.mark, ...))
  }else{
    spntf <- paste0("%.", digits, "f")
    return(sprintf(spntf, round(x, round)))
  }
}

ucid_analysis <- function(input, pct_ID = 0.95){
  
  uniq_MID_grl <- input %$%
    cigarRanges(rname, cigar, pos, sym = "MID", seq.info = ref_seqinfo) %>%
    unlist() %>%
    split(., .$symbol)
  
  if( "I" %in% names(uniq_MID_grl) ){
    
    uniq_I_grl <- split(uniq_MID_grl[["I"]], uniq_MID_grl[["I"]]$index)
    
    uniq_I_idx <- GenomicRanges::findOverlaps(
        norm_targets, uniq_I_grl, maxgap = config$maxDistFromEdit
      ) %>%
      S4Vectors::subjectHits() %>%
      uniq_I_grl[.] %>%
      names(.) %>%
      as.numeric()
    
  }else{
    
    uniq_I_idx <- numeric()
    
  }
  
  if( "D" %in% names(uniq_MID_grl) ){
    
    uniq_D_grl <- split(uniq_MID_grl[["D"]], uniq_MID_grl[["D"]]$index)
    
    uniq_D_idx <- GenomicRanges::findOverlaps(
        norm_targets, uniq_D_grl, maxgap = config$maxDistFromEdit
      ) %>%
      S4Vectors::subjectHits() %>%
      uniq_D_grl[.] %>%
      names(.) %>%
      as.numeric()
    
  }else{
    
    uniq_D_idx <- numeric()
    
  }
  
  if("M" %in% names(uniq_MID_grl)){
    
    uniq_M_grl <- split(uniq_MID_grl[["M"]], uniq_MID_grl[["M"]]$index)
    
    uniq_M_grl <- uniq_M_grl[
      !as.numeric(names(uniq_M_grl)) %in% c(uniq_I_idx, uniq_D_idx)
    ]
    
    widths_M <- sum(width(uniq_M_grl))
    mismatches_M <- str_count(input$md[as.numeric(names(uniq_M_grl))], "[ATGC]")
    identity_M <- (widths_M - mismatches_M) / widths_M
    uniq_M_idx <- as.numeric(names(identity_M))[identity_M >= pct_ID]
    
  }else{
    
    uniq_M_idx <- numeric()
    
  }
  
  input %>%
    dplyr::mutate(
      is.com = seq_len(n()) %in% uniq_M_idx,
      is.ins = seq_len(n()) %in% uniq_I_idx,
      is.del = seq_len(n()) %in% uniq_D_idx,
      is.unc = !is.com & !is.ins & !is.del
    ) %>%
    dplyr::group_by(specimen, rname) %>%
    dplyr::summarise(
      total.cnt = sum(count),
      unc.freq = sum(count[is.unc] / sum(count)),
      com.freq = sum(count[is.com]) / sum(count),
      in.freq = sum(count[is.ins]) / sum(count),
      del.freq = sum(count[is.del] / sum(count))) %>%
    dplyr::ungroup()
  
}

calc_coverage <- function(gr, resolution){
  #Set up coverage gr
  strandless <- gr
  strand(strandless) <- "*"
  gr_ranges <- range(strandless)
  
  window_seqs <- lapply(gr_ranges, function(chr, res){
    seq(start(chr), end(chr), res)
  }, res = resolution)
  
  coverage_grl <- GRangesList(lapply(
    1:length(gr_ranges), function(i, gr_ranges, window_seqs){
      seqname <- seqnames(gr_ranges[i])
      window <- window_seqs[[i]]
      GRanges(
        seqnames = rep(seqname, length(window)),
        ranges = IRanges(
          start = window, width = rep(resolution, length(window))),
        strand = rep("*", length(window)))
    }, gr_ranges = gr_ranges, window_seqs = window_seqs))
  
  coverage_pos <- coverage_grl
  coverage_pos <- GRangesList(lapply(coverage_pos, function(x){
    strand(x) <- rep("+", length(x))
    x}))
  coverage_neg <- coverage_grl
  coverage_neg <- GRangesList(lapply(coverage_pos, function(x){
    strand(x) <- rep("-", length(x))
    x}))
  
  bind_rows(lapply(1:length(coverage_grl), function(i, gr){
    as.data.frame(coverage_grl[[i]], row.names = NULL) %>%
      select(seqnames, start, end, width) %>%
      mutate(
        readCountsPos = countOverlaps(coverage_pos[[i]], gr),
        readCountsNeg = countOverlaps(coverage_neg[[i]], gr)) %>%
      arrange(seqnames)
  }, gr = gr))
}

#' Trim regions with zero edits from the most 5' and 3' ends of the coverage
#' 
#' @usage trim_0_edges(x, tar_pos, min_dist = 50L)
#' 
#' @param x numeric/integer vector representing the coverage at each position
#' from 5' to 3' of the sequence.
#' 
#' @param tar_pos numeric/integer value specifying the position of the targeted
#' site for investigation.
#' 
#' @param min_dist integer value specifying the minimum distance to trim with 
#' respect to the target position, both + and - directions. Example, if target
#' position is at 75, then only regions up to 25 and after 125 will be 
#' considered for trimming (if coverage is 0).

trim_0_edges <- function(x, tar_pos, min_dist){
  
  idx <- seq_along(x)
  is_0 <- S4Vectors::Rle(x == 0)
  not_protected <- idx <= (tar_pos - min_dist) | idx >= (tar_pos + min_dist)
  trim_up <- ifelse(is_0@values[1] == TRUE, is_0@lengths[1], 0)
  
  trim_down <- ifelse(
    is_0@values[length(is_0@values)] == TRUE, 
    is_0@lengths[length(is_0@values)], 
    0
  )
  
  if( length(is_0@values) > 1 ){
    
    trim_regions <- S4Vectors::Rle(
      values = c(TRUE, FALSE, TRUE), 
      lengths = c(trim_up, length(idx) - (trim_up + trim_down), trim_down)
    )
  
    return(as.vector(is_0 & not_protected & trim_regions))
    
  }else{
    
    return(as.vector(is_0 & not_protected))
    
  }
  
}

plot_indel_cov <- function(gr, target, totals, min_dist = 50L){
  
  df <- as.data.frame(GenomicRanges::coverage(gr)) %>% 
    dplyr::mutate(
      specimen = unique(gr$specimen)
    ) %>%
    dplyr::group_by(specimen, group_name) %>%
    dplyr::mutate(
      pos = seq_len(n()),
      trim = !trim_0_edges(
        x = value, 
        tar_pos = start(target)[match(unique(group_name), seqnames(target))],
        min_dist = min_dist
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(trim) %>%
    dplyr::select(-group, -trim) %>%
    dplyr::rename(rname = group_name) %>%
    dplyr::left_join(totals, by = c("specimen", "rname")) %>%
    dplyr::mutate(coverage = (count - value) / count) %>%
    dplyr::filter(rname %in% unique(seqnames(gr)))
  
  target_sites <- as.data.frame(target, row.names = NULL) %>%
    dplyr::filter(seqnames %in% df$rname) %>%
    dplyr::mutate(pos = start) %>%
    dplyr::select(seqnames, pos)
  
  df <- df %>%
    dplyr::mutate(
      rname = factor(rname, levels = target_sites$seqnames)
    ) %>%
    dplyr::arrange(rname) %>%
    dplyr::mutate(
      pos = pos - target_sites$pos[match(rname, target_sites$seqnames)],
      rname = paste0(as.character(rname), " (read depth = ", pNums(count), ")"),
      rname = factor(rname, levels = unique(rname))
    )
  
  df_mins <- df %>%
    dplyr::group_by(specimen, rname) %>%
    dplyr::summarise(
      min_x = min(pos),
      max_x = max(pos),
      min_cov = min(coverage),
      min_cov_format = paste0(
        pNums(100 * (1-min(coverage)), round = 1, digits = 1), 
        "%"
      )
    )

  ggplot() + 
    geom_bar(
      data = df,
      aes(x = pos, y = coverage), 
      fill = "grey40", width = 1L, stat = "identity"
    ) +
    geom_hline(
      data = df_mins, aes(yintercept = min_cov), color = "grey75"
    ) +
    geom_text(
      data = dplyr::filter(df_mins, min_cov < 0.5), 
      aes(x = min_x, y = min_cov, label = min_cov_format),
      hjust = 0, vjust = 0, nudge_x = 2.5, nudge_y = 0.02,
      color = "white", fontface = "bold"
    ) +
    geom_text(
      data = dplyr::filter(df_mins, min_cov >= 0.5), 
      aes(x = min_x, y = min_cov, label = min_cov_format),
      hjust = 0, vjust = 1, nudge_x = 2.5, nudge_y = -0.02,
      color = "white", fontface = "bold"
    ) +
    scale_y_continuous(breaks = c(0, 1), labels = c(1, 0)) +
    facet_wrap(
      ~ rname, ncol = nrow(target), nrow = 1, scales = "free", drop = TRUE
    ) +
    labs(
      x = "Edit Position", 
      y = "Proportion Edited", 
      title = paste0("Specimen: ", unique(gr$specimen))
    )
  
}

split_char <- function(x, len, sep){
  unname(sapply(x, function(y){
    paste(unlist(lapply(
      split(unlist(str_split(y, "")), ceiling(1:nchar(y)/len)), 
      paste, collapse = "")), collapse = sep) }))
}

# Load Supporting information --------------------------------------------------
## Load config file
config <- yaml::yaml.load_file(args$config)

## Sample Info
sample_info <- read.csv(
    file.path(root_dir, config$Sample_Info)
  ) %>%
  dplyr::mutate(specimen = str_extract(sampleName, "[\\w]+"))

if( !is.null(config$Supplemental_Info) ){
  supp_info <- read.csv(
    file.path(root_dir, config$Supplemental_Info)
  )
}

## Sequence Reference
ref_seqs <- Biostrings::readDNAStringSet(args$ref, format = "fasta")

ref_seqinfo <- GenomeInfoDb::Seqinfo(
  seqnames = names(ref_seqs), 
  seqlengths = Biostrings::width(ref_seqs), 
  isCircular = rep(FALSE, length(ref_seqs)), 
  genome = rep(config$RefGenome, length(ref_seqs))
)

base_gen <- getGenome(config$RefGenome)

## Panel targets
panel_targets <- read.csv(
  file.path(root_dir, config$Panel_Path)
)

targets <- structure(panel_targets$locus, names = panel_targets$target)

norm_targets <- normalizeTargets(
  targets, base_gen, ref_seqs, ref_seqinfo, return = "granges"
)

## Get versioning
soft_version <- as.character(read.delim(
  file = file.path(root_dir, ".version"), header = FALSE
))

build_version <- list.files(file.path(root_dir, "etc")) %>%
  grep("build.b", ., value = TRUE) %>%
  stringr::str_extract("b[0-9]+\\.[0-9]+.[0-9]+")

# Load Alignment data ----------------------------------------------------------
message("Loading alignment data.")

input_uniq <- data.table::fread(
    args$unique, header = TRUE, data.table = FALSE) %>%
  dplyr::mutate(
    samplename = stringr::str_extract(qname, "[\\w\\-]+"),
    specimen = stringr::str_extract(samplename, "[\\w]+"),
    edit = panel_targets$edit[match(rname, panel_targets$target)],
    condition = factor(
      sample_info$editing[match(samplename, sample_info$sampleName)], 
      levels = unique(sample_info$editing)
    )
  ) %>%
  dplyr::filter(samplename %in% sample_info$sampleName)

input_chim <- data.table::fread(
    args$chimera, header = TRUE, data.table = FALSE) %>%
  dplyr::mutate(
    samplename = stringr::str_extract(qname, "[\\w\\-]+"),
    specimen = stringr::str_extract(samplename, "[\\w]+"),
    pri.edit = panel_targets$edit[match(pri.seq, panel_targets$target)],
    sec.edit = panel_targets$edit[match(sec.seq, panel_targets$target)],
    condition = factor(
      sample_info$editing[match(samplename, sample_info$sampleName)], 
      levels = unique(sample_info$editing)
    )
  ) %>%
  dplyr::filter(samplename %in% sample_info$sampleName)


# Analysis ---------------------------------------------------------------------
message("Beginning analysis.")

## Specimen Summary -------------------
summary_tbl <- input_uniq %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(
    reps = dplyr::n_distinct(samplename),
    total_reads = sum(count),
    indel_reads = sum(count[tar.indel]),
    on_tar_reads = sum(count[edit == "on"]),
    off_tar_reads = sum(count[edit == "off"]),
    null_reads = sum(count[edit == "null"])
  ) %>%
  dplyr::mutate(
    specimen = factor(
      specimen, 
      levels = unique(str_extract(sample_info$sampleName, "[\\w]+"))
    )
  ) %>%
  dplyr::arrange(specimen)

names(summary_tbl) <- c(
  "Specimen\n ", "Replicates\n ", 
  paste0(
    c("Total", "InDel", "On-target", "Off-target", "Control-Sites"), 
    "\n (Reads)")
)

## On-target Summary ------------------
if( any(input_uniq$edit == "on") ){
  on_tar_summary_tbl <- input_uniq %>%
    dplyr::filter(edit == "on") %>%
    ucid_analysis()
}

## Off-target Summary ------------------
if( any(input_uniq$edit == "off") ){
  off_tar_summary_tbl <- input_uniq %>%
    dplyr::filter(edit == "off") %>%
    ucid_analysis()
}
  
## Null-target Summary ----------------
if( any(input_uniq$edit == "null") ){
  null_tar_summary_tbl <- input_uniq %>%
    dplyr::filter(edit == "null") %>%
    ucid_analysis()
}

## InDel Profiles ---------------------
### Isolate all reads with InDels in the correct regions
mut_uniq_gr <- input_uniq %$%
  cigarRanges(rname, cigar, pos, sym = "ID", seq.info = ref_seqinfo) %>%
  unlist()

edit_gr <- GenomicRanges::findOverlaps(
    norm_targets, mut_uniq_gr, maxgap = config$maxDistFromEdit
  ) %>%
  S4Vectors::subjectHits() %>% 
  mut_uniq_gr[.]

edit_gr$specimen <- input_uniq$specimen[edit_gr$index]
edit_gr$condition <- input_uniq$condition[edit_gr$index]
edit_gr$edit <- input_uniq$edit[edit_gr$index]
edit_gr$count <- input_uniq$count[edit_gr$index]

edit_grl <- split(edit_gr, edit_gr$symbol)

total_specimen_counts <- input_uniq %>% 
  dplyr::group_by(specimen, rname, condition) %>% 
  dplyr::summarise(count = sum(count)) %>%
  dplyr::ungroup()

# Generate report --------------------------------------------------------------
# Normalize file output path
write(c(), file = args$output)
args$output <- normalizePath(args$output)
unlink(args$output)

output_path <- unlist(strsplit(args$output, "/"))
output_dir <- paste(output_path[seq_len(length(output_path)-1)], collapse = "/")
output_file <- output_path[length(output_path)]

if( args$format == "html" & !str_detect(output_file, ".html$") ){
  output_file <- paste0(output_file, ".html")
}

if( args$format == "pdf" & !str_detect(output_file, ".pdf$") ){
  output_file <- paste0(output_file, ".pdf")
}

if( args$data ){
  
  if(args$format == "html"){
    
    save.image(file = file.path(
      output_dir, stringr::str_replace(output_file, ".html$", ".RData")
    ))
    
  }else if( args$format == "pdf" ){
    
    save.image(file = file.path(
      output_dir, stringr::str_replace(output_file, ".pdf$", ".RData")
    ))
    
  }
  
}

figure_path <- file.path(
  output_dir, gsub("[\\w]+$", "figures", output_file, perl = TRUE)
)

null <- dir.create(figure_path)

if( args$tables ){
  
  tables_path <- file.path(
    output_dir, gsub("[\\w]+$", "tables", output_file, perl = TRUE)
  )
  
  null <- dir.create(tables_path)
  
}

if( args$format == "html" ){
  
  rmarkdown::render(
    input = template_path,
    output_format = output_format, 
    output_file = output_file,
    output_dir = output_dir,
    output_options = list(
      "css" = file.path(code_dir, "report_templates/report.css")
    )
  )
  
}else{
  
  rmarkdown::render(
    input = template_path,
    output_format = output_format, 
    output_file = output_file,
    output_dir = output_dir
  )
  
}


if( !args$figures ){
  
  tmp_fig_paths <- c(
    list.files(
      path = figure_path, pattern = "target_editing", full.names = TRUE
    ),
    list.files(
      path = figure_path, pattern = "target_del_profile", full.names = TRUE
    )
  )
  
  cat(sprintf("Removing temporary files: %s\n", tmp_fig_paths), sep = "")
  null <- file.remove(tmp_fig_paths)
  cat("Removing temorary directory:", figure_path, "\n")
  null <- file.remove(figure_path)
  
}

q()
