#!/usr/bin/env nextflow
/*

============================================================================
  BatchBlaster: Nextflow-based BLAST pipeline
============================================================================
  Version: v0.1
  License: Apache-2.0
  Github : https://github.com/vmikk/BatchBlaster
  Author : Vladimir Mikryukov
----------------------------------------------------------------------------
*/


// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Pipeline version
version = '0.1'

// Input and output direcotries
params.input  = false
params.outdir = "${launchDir}/results"

params.method = "blast"                // other methods not implemented yet

// Taxonomy annotation with BLAST
params.blast_taxdb        = false
params.blastdb_fasta      = false
params.blast_task         = "blastn"   // or "megablast"
params.blast_chunksize    = 4000
params.blast_maxts        = 10
params.blast_hsps         = 1
params.blast_wordsize     = false  // 7
params.blast_evalue       = false  // 0.001
params.blast_reward       = false  // 1
params.blast_penalty      = false  // -1
params.blast_gapopen      = false  // 1
params.blast_gapextend    = false  // 2
params.blast_percidentity = false  // 80

if(params.blast_taxdb){
  bastdb_name = file(params.blast_taxdb).name
  bastdb_dir = file(params.blast_taxdb).parent
}

// Check if input path was provided
if (params.input == false) {
  println( "Please provide the input file with sequences in FASTA format with `--input` parameter.")
  exit(1)
}

// Taxonomy annotation
// use globally-dereplicated sequences
process blast {

    label "main_container"

    // publishDir "BLAST_results", mode: 'symlink'
    // cpus 10

    input:
      path input
      path taxdb_dir

    output:
      path "${input.getBaseName()}.m8.gz", emit: blastchunks

    script:

    // Optional parameters
    wordsize  = params.blast_wordsize     ? "-word_size ${params.blast_wordsize}"         : ""
    evalue    = params.blast_evalue       ? "-evalue    ${params.blast_evalue}"           : ""
    reward    = params.blast_reward       ? "-reward    ${params.blast_reward}"           : ""
    penalty   = params.blast_penalty      ? "-penalty   ${params.blast_penalty}"          : ""
    gapopen   = params.blast_gapopen      ? "-gapopen   ${params.blast_gapopen}"          : ""
    gapextend = params.blast_gapextend    ? "-gapextend ${params.blast_gapextend}"        : ""
    gapextend = params.blast_percidentity ? "-perc_identity ${params.blast_percidentity}" : ""
    
    """
    echo -e "Taxonomy annotation with BLAST\n"
    echo -e "Input file: " ${input}
    echo -e "Database  : " ${taxdb_dir}

    blastn \
      -task ${params.blast_task} \
      -query ${input} \
      -db ${taxdb_dir}/${bastdb_name} \
      -strand both \
      -outfmt=6 \
      -max_target_seqs ${params.blast_maxts} \
      -max_hsps ${params.blast_hsps} \
      -out ${input.getBaseName()}.m8 \
      -num_threads ${task.cpus} \
      ${wordsize} \
      ${evalue} \
      ${reward} \
      ${penalty} \
      ${gapopen} \
      ${gapextend} \
      ${gapextend}

    ## Compress results
    gzip -7 ${input.getBaseName()}.m8

    echo "..Done"

    """
}


// Aggregate BLAST results
process blast_merge {

    label "main_container"

    publishDir "BLAST_results", mode: 'symlink'
    // cpus 1

    input:
      path input

    output:
      path "Blast_hits.m8.gz", emit: m8

    script:
    """

    echo -e "Aggregating BLAST hits\n"
    
    cat *.m8.gz > Blast_hits.m8.gz

    echo "..Done"

    """
}

// Workflow
workflow {

  // Input file with sequences (FASTA)
  ch_inp = Channel.fromPath(params.input)

  // Split FASTA sequences into multiple chunks
  ch_fasta = ch_inp.splitFasta(by: params.blast_chunksize, file:true)

  // Run taxonomy annotation
  blast(ch_fasta, bastdb_dir)

}

