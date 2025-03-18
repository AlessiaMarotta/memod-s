# `memod-s`

> **Note**
> A standardised workflow to explore and analyse prokaryotic methylation patterns for Nanopore sequencing data

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) ![C](https://img.shields.io/badge/c-%2300599C.svg?style=for-the-badge&logo=c&logoColor=white) ![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white) ![Bash](https://img.shields.io/badge/bash-%234EAA25.svg?style=for-the-badge&logo=gnu-bash&logoColor=white) ![Snakemake](https://img.shields.io/badge/Snakemake-svg?style=for-the-badge&logo=c&logoColor=white) 


Hi :blush:

`memod-s` is an easy-to-use Snakemake-based workflow that integrates multiple state-of-the-art tools to perform a comprehensive analysis of nanopore sequencing data.
The workflow covers all critical steps, from basecalling, quality control, genome assembly, and annotation to methylation analysis.
It combines genome-wide methylation profiles with genomic features, providing insights into bacterial methylation patterns and their potential biological implications.

## ⚙️ Installation

To set up `memod-s`, you first need to create an environment where Snakemake will be installed. Then, you can clone the repository and run the workflow. 
Create a dedicated environment for `memod-s` using [mamba package manager](https://github.com/mamba-org/mamba):

```
mamba create -n memod-s -c bioconda snakemake
```

This will create an environment called `memod-s` and install Snakemake.

Then, activate the environment:
```
mamba activate memod-s
```
And clone the repo:

```
git clone https://github.com/AlessiaMarotta/memod-s.git
cd memod-s
```
## 🔧 Usage

Run `memod-s` with *--help* or *-h* arguments to see usage instructions:

```
python memod-s --help
```
```
Welcome to memod-s
usage: memod-s [-h] -i INPUT_DIRECTORY [-o OUTPUT_DIRECTORY] [-dm DORADO_MODELS [DORADO_MODELS ...]] [-ml FILTLONG_MIN_LENGTH] [-kp FILTLONG_KEEP_PERCENT]
               [-rr RACON_ROUNDS] [-mr MEDAKA_ROUNDS] [-mm MEDAKA_MODEL] [-ed EGGNOG_DB] [-ab] [-ab_dir ABRICATE_DB_DIR] [-ab_db ABRICATE_DB_NAME] [-q]

Process files for downstream analysis.

options:
  -h, --help            show this help message and exit
  -i INPUT_DIRECTORY, --input_directory INPUT_DIRECTORY
                        Specify the input directory containing fast5 or pod5 files
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Specify the output directory (default: ./memod)
  -dm DORADO_MODELS [DORADO_MODELS ...], --dorado_models DORADO_MODELS [DORADO_MODELS ...]
                        Specify a list of Dorado models
  -ml FILTLONG_MIN_LENGTH, --filtlong_min_length FILTLONG_MIN_LENGTH
                        Specify Filtlong minimum read length
  -kp FILTLONG_KEEP_PERCENT, --filtlong_keep_percent FILTLONG_KEEP_PERCENT
                        Specify Filtlong percent
  -rr RACON_ROUNDS, --racon_rounds RACON_ROUNDS
                        Specify racon rounds
  -mr MEDAKA_ROUNDS, --medaka_rounds MEDAKA_ROUNDS
                        Specify medaka rounds
  -mm MEDAKA_MODEL, --medaka_model MEDAKA_MODEL
                        Specify medaka model
  -ed EGGNOG_DB, --eggnog_db EGGNOG_DB
                        Specify EggNOG database path
  -ab, --abricate       Specify whether to use Abricate
  -ab_dir ABRICATE_DB_DIR, --abricate_db_dir ABRICATE_DB_DIR
                        Specify Abricate database directory
  -ab_db ABRICATE_DB_NAME, --abricate_db_name ABRICATE_DB_NAME
                        Specify Abricate database name
  -q, --quiet           Minimal standard output

```
🛠 Basic usage, to run `memod-s` with default parameters:

```
memod-s -i /path/to/your/input_dir/with/fast5/or/pod5
```
🔬 Advanced Usage, for a full analysis with specific parameters:

```
memod-s -i /path/to/your/input_dir/with/fast5/or/pod5 \
        -o /path/to/output_dir \
        -ml 1500 \
        -kp 80 \
        -rr 5 \
        -mr 2 \
        -mm r1041_e82_400bps_sup_v4.3.0 \
        -dm dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1 \
            dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1 \
        -ab \
        -ab_db vfdb \
        -q
```

## 🐍 Workflow

### Core

1. Convert fast5 files to pod5 with [pod5 convert fast5](https://pod5-file-format.readthedocs.io/en/latest/docs/tools.html#pod5-convert-fast5)
2. Download specific methylation models with Dorado. Check the [avaiable basecalling models](https://github.com/nanoporetech/dorado?tab=readme-ov-file#available-basecalling-models)
3. Basecalling with [Dorado](https://github.com/nanoporetech/dorado)
4. Map reads to reference with [minimap2](https://github.com/lh3/minimap2)
5. Quality check with [NanoPlot](https://github.com/wdecoster/NanoPlot)
6. Quality filtering with [Filtlong](https://github.com/rrwick/Filtlong)
7. Quality check after fitlering with [NanoPlot](https://github.com/wdecoster/NanoPlot)
8. Assembly with [Dragonflye](https://github.com/rpetit3/dragonflye)
9. Assembly evaluation with [quast](https://github.com/ablab/quast)
10. Genome annotation with [prokka](https://github.com/tseemann/prokka)
11. Functional annotation with [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)
12. Call methylation with [MicrobeMod](https://github.com/cultivarium/MicrobeMod)
13. Combine genome-wide methylation profiles with genomic features with [MeStudio](https://github.com/combogenomics/MeStudio) 
14. Circular density plot for each motif with [Circlize](https://github.com/jokergoo/circlize)
15. Gene Set Enrichment Analysis with [fgsea](https://github.com/alserglab/fgsea)

### Bonus

1.  Antimicrobial resistance or virulence genes screening with [Abricate](https://github.com/tseemann/abricate)

## 🗂️ Test memod-s

* You can use the sample data available at this [link] to test the workflow.


## 📄 Publications

The `memod-s` workflow 

```
NOSTRO ARTICOLO
```
## 🖇️ References
```
Riccardi, Christopher, Iacopo Passeri, Lisa Cangioli, Camilla Fagorzi, Marco Fondi, and Alessio Mengoni. 2023.
"Crossing Bacterial Genomic Features and Methylation Patterns with MeStudio: An Epigenomic Analysis Tool"
International Journal of Molecular Sciences 24, no. 1: 159.
https://doi.org/10.3390/ijms24010159

Crits-Christoph, Alexander, Shinyoung Clair Kang, Henry H. Lee, and Nili Ostrov.
"MicrobeMod: A computational toolkit for identifying prokaryotic methylation and restriction-modification with nanopore sequencing."
bioRxiv (2023): 2023-11. doi: https://doi.org/10.1101/2023.11.13.566931
```
