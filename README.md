# BatchBlaster  <img src='assets/BatchBlaster_Logo.webp' align="right" height="200" />

![Version](https://img.shields.io/badge/version-0.0.1-blue)
![License](https://img.shields.io/github/license/vmikk/BatchBlaster)

BatchBlaster is a bioinformatics pipeline that employs BLAST (Basic Local Alignment Search Tool), an essential algorithm for comparing primary biological sequence information, to perform efficient and high-throughput taxonomic identification searches.  

BatchBlaster is built using the Nextflow workflow management system, ensuring portability and reproducibility across multiple platforms.  

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
## Output

The results will be saved in the specified output directory (`./results` by default). It includes:

- BLAST search results in tabular format ([`m8` a.k.a. `-outfmt 6`](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6))  
- Summary report  


## License

This project is licensed under the terms of the Apache-2.0 license.

---

Please feel free to submit issues and pull requests, your contributions are welcome!

