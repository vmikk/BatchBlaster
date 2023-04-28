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

############################################## Load packages and data

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("plyr")
load_pckg("metagMisc")
load_pckg("Biostrings")
load_pckg("openxlsx")

cat("\n")

## Set CPU thread number
cat("Setting number of CPU threads to: ", CPUTHREADS, "\n")
setDTthreads(threads = CPUTHREADS)  # for data.table


##### Parse BLAST

## Load sequences
cat("Loading FASTA sequneces\n")
seqs <- FASTA
seqs <- readDNAStringSet(seqs, format="fasta")
seqs_names <- names(seqs)
names(seqs) <- gsub(pattern = ";.*", replacement = "", x = names(seqs))

## Load reference sequences
cat("Loading BLAST database\n")
db_refs <- DATABASE
db_refs <- readDNAStringSet(db_refs, format="fasta")
## BLAST columns
bcolz <- c(
  "QueryName", "TargetName",
  "SeqIdentity", "AlignLen", "MismatchN", "GapOpenings",
  "QueryStart", "QueryEnd", "TargetStart", "TargetEnd", "Evalue", "BitScore")

## Load blast results
cat("Reading m8 file\n")
BLASTS_10h <- read_m8(INPUT, blast_colz = bcolz, package = "data.table")

## Remove size annotations
cat("Extracting query IDs\n")
BLASTS_10h[ , QueryName := tstrsplit(QueryName, ";", keep=1) ]
  cat("Splitting taxonomy string\n")
  BLASTS_10h[, c("AccID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Function", "LT_comment", "SH") := tstrsplit(TargetName, ";", keep=1:11) ]

  tcolz <- c("AccID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Function", "LT_comment", "SH")
  bcolz <- c(bcolz[-2], tcolz)
  BLASTS_10h <- BLASTS_10h[, ..bcolz]

  ## Convert to wide format
  cat("Converting to wide format\n")
  BW <- blast_to_wide(BLASTS_10h, max_hits = 10,
    taxonomy = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Function", "LT_comment", "SH"),
    seqs = seqs, refs = db_refs)


## Add missing sequences (no BLAST info)
if(any(!names(seqs) %in% BW$QueryName)){
  cat("WARNING: some sequences do not have BLAST annotations\n")

  ## Add missing sequences
  SEQS_TAB <- data.table(
    QueryName = names(seqs),
    Seq = as.character(seqs))

  SEQS_TAB <- SEQS_TAB[ ! QueryName %in% BW$QueryName ]

  cat("..", nrow(SEQS_TAB), " sequences without BLAST hit were added to the wide table\n")
  BW <- rbind(BW, SEQS_TAB, fill = TRUE)
}

