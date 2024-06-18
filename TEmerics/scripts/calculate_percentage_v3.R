#!/usr/bin/env Rscript

#################
###  IMPORTS  ###
#################
# Import required packages
if (
  suppressWarnings(suppressPackageStartupMessages(require("data.table"))) == FALSE
) {
  stop("[ERROR] Package 'data.table' required! Aborted.")
}
if (
  suppressWarnings(suppressPackageStartupMessages(require("dplyr"))) == FALSE
) {
  stop("[ERROR] Package 'dplyr' required! Aborted.")
}
if (
  suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE
) {
  stop("[ERROR] Package 'optparse' required! Aborted.")
}
if (
  suppressWarnings(suppressPackageStartupMessages(require("rstatix"))) == FALSE
) {
  stop("[ERROR] Package 'rstatix' required! Aborted.")
}


#######################
###  PARSE OPTIONS  ###
#######################


# Get script name
script <- sub("--file=", "", basename(commandArgs(trailingOnly = FALSE)[4]))

# Build description message
description <- "ADD DESCRIPTION."
author <- "Author: Marta Coronado <mail marta>" # Can be removed
maintainer <- "Maintainer: Albert Garcia <your mail>" # Can be removed
version <- "Version: 1.1.0 (MAR-2024)"
requirements <- "Requires: data.table, dplyr, optparse and rstatix"
msg <- paste(description, author, maintainer, version, requirements, sep = "\n")

# Define list of arguments
option_list <- list(
  make_option(
    "--in-file",
    action = "store",
    type = "character",
    help = "ADD DESCRIPTION. Required!",
    metavar = "directory"
  ),
  make_option(
    "--ref-trans",
    action = "store",
    type = "character",
    help = "ADD DESCRIPTION. Required!",
    metavar = "directory"
  ),
  make_option(
    "--out-file",
    action = "store",
    type = "character",
    default = file.path(getwd(), "chimerics.filterLengthMin100.filterScovMin80.blastn"),
    help = "ADD DESCRIPTION. Default: $PWD/chimerics.filterLengthMin100.filterScovMin80.blastn",
    metavar = "directory"
  ),
  make_option(
    c("-h", "--help"),
    action = "store_true",
    default = FALSE,
    help = "Show this information and die."
  ),
  make_option(
    c("-u", "--usage"),
    action = "store_true",
    default = FALSE,
    dest = "help",
    help = "Show this information and die."
  ),
  make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Print log messages to STDOUT."
  )
)

# Parse command-line arguments
opt_parser <-
  OptionParser(
    usage = paste(
      "Usage:",
      script,
      "[OPTIONS] --in-file <path/to/input/file>",
      "--ref-trans <path/to/input/reference_transcripts.tab>",
      "--out-file <path/to/output/file>\n",
      sep = " "
    ),
    option_list = option_list,
    add_help_option = FALSE,
    description = msg
  )
opt <- parse_args(opt_parser)

# Re-assign variables
in_file <- opt$`in-file`
ref_trans <- opt$`ref-trans`
out_file <- opt$`out-file`
verb <- opt$`verbose`

# Validate required arguments
if (is.null(in_file) || is.null(ref_trans)) {
  print_help(opt_parser)
  stop("[ERROR] Required argument missing! Aborted.")
}

######################
###      MAIN      ###
######################

# Write log
if (verb) {
  cat(paste("Reading input: ", in_file, "...\n", sep = ""), sep= "")
}

blastn_result <- read.table(in_file)
colnames(blastn_result) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")


# Write log
if (verb) {
  cat("Filtering by length and coverage...\n", sep= "")
}

blastn_result <- blastn_result[blastn_result$length > 100,]
# should work
blastn_result$scov <- ((blastn_result$send - blastn_result$sstart + 1) / blastn_result$length)
blastn_result <- blastn_result[blastn_result$scov>0.8,]


# Write log
if (verb) {
  cat("Filtering by reference transcript and bitscore...\n", sep= "")
}

reference_transcripts <- read.table(ref_trans, header = T)

blastn_result_merge <- merge(blastn_result, reference_transcripts[,c("ref_gene_id", "qry_id")], by.x="sseqid", by.y="qry_id")

blastn_result_merge_filter <- blastn_result_merge %>% group_by(qseqid,ref_gene_id) %>%  filter(bitscore == max(bitscore))
blastn_result_merge_filter<-blastn_result_merge_filter[,c("qseqid","sseqid","ref_gene_id","pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore","scov")]

# Write log
if (verb) {
  cat("Filtering by unique combinations of ID and reference gene ID...\n", sep= "")
}

# .keep_all = TRUE to keep all columns
blastn_result_unique <- blastn_result_merge_filter %>% distinct(qseqid, ref_gene_id, .keep_all = TRUE)

# Write log
if (verb) {
  cat(paste("Writing table: ", out_file, "\n", sep = ""), sep = "")
}

write.table(blastn_result_unique,file=out_file, quote = F, sep = "\t", col.names = F, row.names = F)


# Write log
if (verb) {
  cat("Done.\n", sep = "")
}
