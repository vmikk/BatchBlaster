# BatchBlaster  <img src='assets/BatchBlaster_Logo.webp' align="right" height="200" />

![Version](https://img.shields.io/badge/version-0.0.1-blue)
![Lifecycle](https://img.shields.io/badge/status-experimental-orange)
![License](https://img.shields.io/github/license/vmikk/BatchBlaster)

BatchBlaster is a bioinformatics pipeline that employs BLAST (Basic Local Alignment Search Tool), an essential algorithm for comparing primary biological sequence information, to perform efficient and high-throughput taxonomic identification searches.  

BatchBlaster is built using the [Nextflow](https://www.nextflow.io/) workflow management system, ensuring portability and reproducibility across multiple platforms. The pipeline is primarily designed for use on High Performance Computing (HPC) clusters, including the capability to submit tasks to the SLURM job scheduling system.  

The name 'BatchBlaster' originates from its robust capability to submit and process BLAST tasks in batches, optimizing for speed and performance in large-scale sequence analysis tasks.  

## Features

- High throughput BLAST search  
- Scalable and reproducible analysis with Nextflow  
- Multi-platform compatibility (Linux, MacOS, Windows)  

## Quick Start

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)

    ```bash
    curl -s https://get.nextflow.io | bash
    ```
2. Run BatchBlaster

    ```bash
    nextflow run vmikk/BatchBlaster -r main --input 'path/to/your/input' ...
    ```

## Parameters

- `--input` : Path to the input file containing the sequences (Required)  
- `--outdir` : Path to the output directory (Default: `./results`)  
- `--blast_taxdb` : Path to the BLAST database  
- ...

## Output

The results will be saved in the specified output directory (`./results`, by default). Output includes:

- BLAST search results in tabular format ([`m8` a.k.a. `-outfmt 6`](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6))  
- A table with best BLAST hits reshaped into wide format  
- Summary report  

## Dependencies

- [Nextflow](https://www.nextflow.io/) (>=23.04.0)  
- [Singularity](https://sylabs.io/singularity/) or [Docker](https://www.docker.com/)


## Future Plans

- **Integration of additional sequence analysis methods** (e.g., MMSeqs2, SINTAX, etc.)  
- **Inclusion of Lowest Common Ancestor (LCA) estimation**  
- Implementation of domain-specific threshold filtering for taxonomic annotation (e.g., for fungal sequences)  
- Adding advanced machine learning algorithms for more accurate taxonomic classification (e.g., deep learning models that have been trained on the [UNITE database](https://unite.ut.ee/index.php))  
- Implementation of a hybrid annotation approach (e.g., integration of classification results from various methods to enhance accuracy and reliability of taxonomic identification)  

We are excited to share these enhancements in our forthcoming updates, so stay tuned!

## License

This project is licensed under the terms of the Apache-2.0 license.

---

Please feel free to submit [issues](https://github.com/vmikk/BatchBlaster/issues) 
and [pull requests](https://github.com/vmikk/BatchBlaster/pulls), your contributions are welcome!

