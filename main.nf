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
params.input = false
params.outdir = "${launchDir}/results"

// Taxonomy annotation with BLAST
params.blast_taxdb = false
params.blast_task = "blastn"   // or "megablast"
params.blast_chunksize = 4000
params.blast_maxts = 10
params.blast_hsps = 1
params.blast_wordsize     = false  // 7
params.blast_evalue       = false  // 0.001
params.blast_reward       = false  // 1
params.blast_penalty      = false  // -1
params.blast_gapopen      = false  // 1
params.blast_gapextend    = false  // 2
params.blast_percidentity = false  // 80
