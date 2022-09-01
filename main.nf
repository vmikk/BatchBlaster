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
