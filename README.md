# ONT-DMS

ONT-DMS is a Nextflow pipeline for processing barcoded protein variant libraries and [NestLink libraries](https://www.nature.com/articles/s41592-019-0389-8) sequenced by nanopore sequencing.
Reads are binned according to their barcodes or flycodes (UMIs).
Accurate consensus sequences are calculated using Dorado polish.
Finally, variants are called with the pipeline, linking barcodes or flycodes with their respective protein variants.

> [!NOTE]
> Tested on macOS 26.7 (Mamba), Windows 11 (WSL with Docker), and Ubuntu 24.04 LTS (Mamba, Apptainer).

## Requirements
### Local execution

- Nextflow ([Installation guide](https://www.nextflow.io/docs/latest/install.html)).
- Mamba/ Conda ([https://conda-forge.org/](https://conda-forge.org/)) or Docker
- Dorado ([Installation guide](https://software-docs.nanoporetech.com/dorado/latest/#installation)), does not have to be installed when using containers.

### Cluster execution

- Slurm workload manager
- Nextflow, installed in a Mamba/ Conda environment called `nextflow`. 
- Mamba/ Conda or Apptainer
- Dorado if not using containers

## Running the pipeline

1. Clone the repository with `git clone https://github.com/fabianackle/ONT-DMS.git`.
2. Edit the nextflow configuration file `nextflow.config`.
3. Create a `params.json` file with the parameters listed below, specify the nanopore reads (BAM) and reference sequence, see the examples contained in this repo.
4. Run the pipeline with either `./run_ONT-DMS.sh` for local execution or on a cluster with `sbatch run_ONT-DMS.slurm`.

> [!TIP]
> The pipeline can also be run directly: `nextflow run main.nf -profile conda -params-file params.json`.

# Parameters

| Parameter                 | Type                 | Description                                                                          |
|---------------------------|----------------------|--------------------------------------------------------------------------------------|
| `data`                    | String               | Path to input BAM file(s).†                                                          |
| `reference`               | String               | Path to reference FASTA file.                                                        |
| `filter_quality`          | Float                | Minimum mean read quality threshold.                                                 |
| `filter_min_length`       | Integer              | Read filtering minimum length threshold.                                             |
| `filter_max_length`       | Integer              | Read filtering maximum length threshold.                                             |
| `extract_barcode_adapter` | String               | Linked cutadapt adapter for barcode extraction.†2                                    |
| `barcode_regex`           | String               | Regular expression matching the barcode.                                             |
| `barcode_min_coverage`    | Integer              | Minimum number of reads containing a barcode required to classify it as high quality.|
| `polish_bacteria`         | Boolean              | Use a bacterial model for polishing with dorado (optional).†3                        |
| `barcode_5p`              | String               | 5' sequence flanking the barcode.                                                    |
| `barcode_3p`              | String               | 3' sequence flanking the barcode.                                                    |
| `orf_5p`                  | String               | 5' sequence flanking the ORF.                                                        |
| `orf_3p`                  | String               | 3' sequence flanking the ORF.                                                        |
| `translate_barcode`       | Boolean              | Translates barcode, used with flycodes.                                              |
| `outdir`                  | String               | Output directory for results.                                                        |

† for multiple BAM files use `*`, e.g. `data/barcode*.bam`.
†2 see [Linked adapters (combined 5’ and 3’ adapter)](https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapters).
†3 do not use if you want to use move tables, as currently there are no models that support both.
