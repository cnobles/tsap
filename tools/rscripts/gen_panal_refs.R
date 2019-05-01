# Rscript to generate referece sequences for panel targets.
# Requires: Biostrings, Biogenerics, BSgenome, GenomicRanges, parallel

# Input parsing and setup ------------------------------------------------------
options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))

args <- commandArgs(trailingOnly = FALSE)

code_dir <- dirname(
  sub("--file=", "", grep("--file=", args, value = TRUE)))

#' Set up and gather command line arguments
parser <- ArgumentParser(
  description = "R-based utility to generate reference sequences for target panel.")
parser$add_argument(
  "panelFile", nargs = 1, type = "character",
  help = paste(
    "CSV file with columns: target,locus,fwd.primer,rev.primer.", 
    "Primer sequences are just portion annealing to genomic DNA and", 
    "will be included in sequence data."))
parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", 
  help = "Output file name. Output type will be parsed from extension.")
parser$add_argument(
  "-r", "--ref", nargs = 1, type = "character", default = "hg38",
  help = "Genomic reference draft. Default: hg38.")
parser$add_argument(
  "-m", "--maxLength", nargs = 1, type = "integer", default = 400,
  help = "Maximum length for generated sequences. Default at 400 nts.")
parser$add_argument(
  "-c", "--cores", nargs = 1, type = "integer", default = 1,
  help = "If greater than 1, script will use r-parallel to multithread process.")
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

input_table <- data.frame(
  "Variable" = names(args), 
  "Value" = sapply(
    seq_along(args), function(i) paste(args[[i]], collapse = ", ")))
input_table <- input_table[
  match(
    c("panelFile", "output", "ref", "maxLength", "cores"), 
    input_table$Variable),]
input_table$Variable <- paste0(format(input_table$Variable), " :")
cat("Demultiplex Inputs\n")
print.data.frame(
  data.frame(input_table, row.names = NULL),
  row.names = FALSE, right = FALSE)

# Function definitions ---------------------------------------------------------
seq_file_type <- function(file){
  seqType <- unlist(strsplit(file, "/"))
  seqType <- seqType[length(seqType)]
  seqType <- stringr::str_extract(seqType, ".fa[\\w]*")
  if(any(!seqType %in% c(".fa", ".fq", ".fasta", ".fastq"))){
    stop(paste(
      "Unrecognized sequence file type, please convert to '*.fasta' or", 
      "'*.fastq'. Gzip compression is acceptable as well."))
  }else{
    return(ifelse(seqType %in% c(".fa", ".fasta"), "fasta", "fastq"))
  }
}

get_genome <- function(ref){
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
  
isolate_region <- function(fwd.primer, rev.primer, target, ref, 
                           max.size = 400, drop.alt.chr = TRUE){
  if(drop.alt.chr){
    acceptable_seqs <- GenomeInfoDb::seqlevels(ref)[
      !stringr::str_detect(GenomeInfoDb::seqlevels(ref), "_")]
  }else{
    acceptable_seqs <- GenomeInfoDb::seqlevels(ref)
  }
  
  fwd_reg <- Biostrings::vmatchPattern(fwd.primer, ref)
  rev_reg <- Biostrings::vmatchPattern(rev.primer, ref)
  
  fwd_reg <- fwd_reg[GenomicRanges::seqnames(fwd_reg) %in% acceptable_seqs]
  rev_reg <- rev_reg[GenomicRanges::seqnames(rev_reg) %in% acceptable_seqs]
  
  fwd_reg$type <- "fwd"
  rev_reg$type <- "rev"
  
  # Pair regions together and check for approriate pairing
  regions <- c(fwd_reg, rev_reg)
  red_regions <- GenomicRanges::reduce(
    regions, min.gapwidth = max.size, ignore.strand = TRUE, with.revmap = TRUE)
  red_regions <- red_regions[sapply(
    red_regions$revmap, function(x) all(c("fwd", "rev") %in% regions$type[x]))]
  red_regions <- red_regions[width(red_regions) <= max.size]
  red_regions <- red_regions[
    GenomicRanges::seqnames(red_regions) == 
        stringr::str_extract(target, "[\\w]+") &
      BiocGenerics::start(red_regions) < 
        as.numeric(stringr::str_extract(target, "[0-9]+$")) &
      BiocGenerics::end(red_regions) > 
        as.numeric(stringr::str_extract(target, "[0-9]+$")) ]
  
  # Return sequence associated with region
  unique(BSgenome::getSeq(ref, red_regions))
}

# Load data and isolate priming sequences --------------------------------------
if(seq_file_type(args$output) == "fastq"){
  stop("Please change output file name to fasta format.")
}

message("\nLoading reference genome: ", args$ref)

ref_genome <- get_genome(args$ref)

if(exists("ref_genome")) message("\nReference genome loaded.")

panel_data <- read.csv(args$panelFile)

cat("\nTarget Panel Table\n")
print.data.frame(panel_data, row.names = FALSE, width = Inf)

message("\nIsolating targeted region sequences.")

if(args$cores <= 1){
  panel_seqs <- structure(lapply(
    seq_len(nrow(panel_data)), function(i, panel, ref, args){
      isolate_region(
        fwd.primer = panel$fwd.primer[i], 
        rev.primer = panel$rev.primer[i], 
        target = panel$locus[i], 
        ref = ref,
        max.size = args$maxLength)
    }, panel = panel_data, ref = ref_genome, args = args), 
    names = panel_data$target)
}else{
  buster <- parallel::makeCluster(args$cores)
  
  parallel::clusterExport(buster, varlist = c("args", "isolate_region"))
  
  panel_seqs <- structure(lapply(
    seq_len(nrow(panel_data)), function(i, panel, ref, args){
    isolate_region(
      fwd.primer = panel$fwd.primer[i], 
      rev.primer = panel$rev.primer[i], 
      target = panel$locus[i], 
      ref = ref,
      max.size = args$maxLength)
    }, panel = panel_data, ref = ref_genome, args = args), 
    names = panel_data$target)

  parallel::stopCluster(buster)
}

if(exists("panel_seqs")) message("Sequences for targeted regions obtained.")

# Write output -----------------------------------------------------------------
output_seqs <- unlist(DNAStringSetList(panel_seqs))

# Write sequences to output file.
Biostrings::writeXStringSet(
  x = output_seqs, 
  filepath = args$output, 
  append = FALSE, 
  format = "fasta",
  width = max(c(width(output_seqs), 1)))

if(file.exists(args$output)){
  message("Output writen to file: ", args$output)
}else{
  message("Unable to verify output file's existance.")
}
message("\nCompleted script.")
q()