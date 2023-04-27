#!/usr/bin/env Rscript

## Script to parse BLAST results

############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-m", "--m8"),    action="store", default=NA, type='character', help="BLAST results in m8 format"),
  make_option(c("-f", "--fasta"), action="store", default=NA, type='character', help="Sequences that were BLASTed (FASTA)"),
  make_option(c("-d", "--db"),    action="store", default=NA, type='character', help="BLAST database (FASTA)"),
  
  make_option(c("-u", "--outputprefix"), action="store", default="Blast_hits", type='character', help="Output file prefix"),
  make_option(c("-t", "--threads"),      action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4")
)
opt <- parse_args(OptionParser(option_list=option_list))

## Assign variables
INPUT      <- opt$m8
FASTA      <- opt$fasta
DATABASE   <- opt$db
OUTPUT     <- opt$outputprefix
CPUTHREADS <- as.numeric( opt$threads )

## Log assigned variables
cat(paste("BLAST results (m8): ", INPUT,    "\n", sep=""))
cat(paste("FASTA sequences: ",    FASTA,    "\n", sep=""))
cat(paste("Database: ",           DATABASE, "\n", sep=""))
cat(paste("Output prefix: ",                 OUTPUT,     "\n", sep=""))
cat(paste("Number of CPU threads to use: ",  CPUTHREADS, "\n", sep=""))

cat("\n")

