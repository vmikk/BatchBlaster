#!/usr/bin/env Rscript

## Script to parse BLAST results

## Usage:
# ./parse_blast_results.R \
#    --m8           Blast_hits.m8.gz \
#    --fasta        tst.fa \
#    --db           UNITE_9.7.fasta.gz \
#    --maxhits      10 \
#    --outputprefix Blast_hits \
#    --threads      4 \
#    --splittax     TRUE \
#    --taxcolumns   "AccID,Kingdom,Phylum,Class,Order,Family,Genus,Species,Function,Dataset" \
#    --exportexcel  TRUE

## !NB. If `--splittax == TRUE`,
##      `--taxcolumns` should be a comma-separated list
##      with the first element "AccID" (accession ID)
##      In the taxonomy database, headers should have SEMICOLON (;) as field separator!
## 
## !NB. Export in Excel format works only if the table has < 200 000 rows


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
  
  make_option(c("-w", "--maxhits"),    action="store", default=10, type='integer', help="Maximum number of BLAST hits to keep (default, 10)"),
  make_option(c("-s", "--splittax"),   action="store", default=TRUE, type='logical', help="Logical, split taxonomy string (default, TRUE"),
  make_option(c("-x", "--taxcolumns"), action="store", default="AccID,Kingdom,Phylum,Class,Order,Family,Genus,Species,Function,Dataset", type='character', help="Field names of the database (comma-separated)"),

  make_option(c("-u", "--outputprefix"), action="store", default="Blast_hits", type='character', help="Output file prefix"),
  make_option(c("-e", "--exportexcel"),  action="store", default=TRUE, type='logical', help="Export results in XLSX format (could fail for large number of records)"), 
  make_option(c("-t", "--threads"),      action="store", default=4L, type='integer', help="Number of CPU threads for arrow, default 4")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Validation of the required argiments
if(is.na(opt$m8)){
  cat("m8 is not specified.\n", file=stderr())
  stop()
}
if(is.na(opt$fasta)){
  cat("Input sequences are not provided.\n", file=stderr())
  stop()
}
if(is.na(opt$db)){
  cat("Sequence database not specified.\n", file=stderr())
  stop()
}

## Assign variables
INPUT      <- opt$m8
FASTA      <- opt$fasta
DATABASE   <- opt$db
MAXHITS    <- as.integer( opt$maxhits )
SPLITTAX   <- as.logical( opt$splittax )
TAXCOLS    <- opt$taxcolumns
OUTPUT     <- opt$outputprefix
EXPORTXLSX <- as.logical( opt$exportexcel )
CPUTHREADS <- as.numeric( opt$threads )

## Log assigned variables
cat(paste("BLAST results (m8): ", INPUT,    "\n", sep=""))
cat(paste("FASTA sequences: ",    FASTA,    "\n", sep=""))
cat(paste("Database: ",           DATABASE, "\n", sep=""))
cat(paste("Max hits to keep: ",   MAXHITS,  "\n", sep=""))
cat(paste("Split tax info: ",     SPLITTAX, "\n", sep=""))
cat(paste("Tax database column names: ",     TAXCOLS,    "\n", sep=""))
cat(paste("Output prefix: ",                 OUTPUT,     "\n", sep=""))
cat(paste("Export Excel (xlsx) files: ",     EXPORTXLSX, "\n", sep=""))

cat(paste("Number of CPU threads to use: ",  CPUTHREADS, "\n", sep=""))

cat("\n")



############################################## Data for debuging

# INPUT      <- "Blast_hits.m8.gz"
# FASTA      <- "inp.fasta"
# DATABASE   <- "UNITE_9.7_230310.fasta.gz"
# OUTPUT     <- "Blast_hits"
# EXPORTXLSX <- TRUE
# CPUTHREADS <- 4
# MAXHITS    <- 10
# SPLITTAX   <- TRUE
# TAXCOLS    <- "AccID,Kingdom,Phylum,Class,Order,Family,Genus,Species,Function,Dataset"


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
if(EXPORTXLSX == TRUE){ load_pckg("openxlsx") }

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
cat("..Number of sequences detected: ", length(seqs), "\n")

## Check if input sequence names are unique
if(length(names(seqs)) != length(unique(names(seqs)))){
  cat("..WARNING: Input sequence names are not unique!\n")
  cat("..Therefore, merging with taxonomic information could be incorrect.\n")
}

## Load reference sequences
cat("Loading BLAST database\n")
db_refs <- DATABASE
db_refs <- readDNAStringSet(db_refs, format="fasta")
cat("..Number of records in the database: ", length(db_refs), "\n")

if(SPLITTAX == TRUE){
  ## Keep only Accession ID for sequences in the databse
  names(db_refs) <- do.call(rbind, strsplit(x = names(db_refs), split = ";"))[,1]
}

## BLAST columns
bcolz <- c(
  "QueryName", "TargetName",
  "SeqIdentity", "AlignLen", "MismatchN", "GapOpenings",
  "QueryStart", "QueryEnd", "TargetStart", "TargetEnd", "Evalue", "BitScore")

## Load blast results
cat("Reading m8 file\n")
BLASTS_10h <- read_m8(INPUT, blast_colz = bcolz, package = "data.table")

cat("..Number of BLAST hits in total: ", nrow(BLASTS_10h), "\n")

## Remove size annotations
cat("Extracting query IDs\n")
BLASTS_10h[ , QueryName := tstrsplit(QueryName, ";", keep=1) ]

## Split tax by tax ranks
if(SPLITTAX == TRUE){

  cat("Splitting taxonomy string\n")

  ## Prepare a vector of column names from the database
  tcolz <- strsplit(x = TAXCOLS, split = ",")[[1]]
  cat(".. ", length(tcolz), " column names are provided for the database\n")  

  ## Split headers of the matches
  ## NB. extra columns (not specified with TAXCOLS) will be ignored
  cat(".. Splitting databse matches\n")
  BLASTS_10h[, eval(tcolz) := tstrsplit(TargetName, ";", keep=1:length(tcolz)) ]

  ## Remove non-splited columns ("QueryName", "TargetName") 
  bcolz <- c(bcolz[-2], tcolz)
  BLASTS_10h <- BLASTS_10h[, ..bcolz]

  ## Convert to wide format
  cat("Converting to wide format\n")
  BW <- blast_to_wide(
    BLASTS_10h,
    max_hits = MAXHITS,
    taxonomy = tcolz[-1],          # all columns except `AccID`
    seqs = seqs, refs = db_refs)

} else {
  ## Do not split taxonomy, keep target ID as-is

  cat("Converting to wide format\n")
  BW <- blast_to_wide(BLASTS_10h,
    max_hits = MAXHITS,
    taxonomy = "TargetName",
    seqs = seqs, refs = db_refs)

}



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


# Remove redundant seqs
if(length(unique(BW$QueryName)) != nrow(BW)){
  cat("WARNING: there are some duplicated query sequences (based on ID); they will be removed\n")
  BW <- unique(BW, by = "QueryName")
}



cat("Exporting RData\n")

## Test if pigz is available
pigz <- Sys.which("pigz")

if(! pigz %in% ""){

  cat("..Exporting with pigz [parallel gzip]")

  ## Parallel save
  saveRDS.gz <- function(object, file, threads = 2) {
    con <- pipe(paste0("pigz -p", threads," > ",file),"wb")
    saveRDS(object, file = con)
    close(con)
  }

  cat("..Wide table\n")
  saveRDS.gz(
    object  = BW,
    file    = paste0(OUTPUT, "_wide.RData"),
    threads = CPUTHREADS)

  cat("..Long table\n")
  saveRDS.gz(
    object  = BLASTS_10h,
    file    = paste0(OUTPUT, "_long.RData"),
    threads = CPUTHREADS)


} else {

  cat("..Wide table\n")
  saveRDS(
    object   = BW,
    file     = paste0(OUTPUT, "_wide.RData"),
    compress = "xz")

  cat("..Long table\n")
  saveRDS(
    object   = BLASTS_10h,
    file     = paste0(OUTPUT, "_long.RData"),
    compress = "xz")

} # end of RData export



## Export results in Excel format
if(EXPORTXLSX == TRUE){ 
  cat("Exporting Excel table\n")
  cat("..There are ", nrow(BW), "rows in the table\n")

  if(nrow(BW) < 200000){
  
    write.xlsx(list(
      "BLAST" = BW
      ),
      file = paste0(OUTPUT, "_BestHits.xlsx"),
      colNames = TRUE)
  
  } else {

    cat("..The table is probably too large for Excel\n")
    cat("..Therefore, writing data in tab-separated text file\n")

    fwrite(x = BW,
      file = paste0(OUTPUT, "_BestHits.txt.gz"),
      sep = "\t", compress = "gzip")

  }
} # end of tabular export


cat("All done!\n")

##################### Session info

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")

cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
