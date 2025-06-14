# 🧬 `memod-s`

> **Note**
> A standardised workflow to explore and analyse prokaryotic methylation patterns for Nanopore sequencing data

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) ![C](https://img.shields.io/badge/c-%2300599C.svg?style=for-the-badge&logo=c&logoColor=white) ![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white) ![Bash](https://img.shields.io/badge/bash-%234EAA25.svg?style=for-the-badge&logo=gnu-bash&logoColor=white) ![Snakemake](https://img.shields.io/badge/Snakemake-svg?style=for-the-badge&logo=c&logoColor=white) 



`memod-s` is an easy-to-use Snakemake-based workflow that integrates multiple state-of-the-art tools to perform a comprehensive analysis of nanopore sequencing data.
The workflow covers all critical steps, from basecalling, quality control, genome assembly, and annotation to methylation analysis.
It combines genome-wide methylation profiles with genomic features, providing insights into bacterial methylation patterns and their potential biological implications.

![memod-s scheme](https://github.com/user-attachments/assets/49da3727-1f49-4349-bdf2-d44dfde4aba4)

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
./memod-s --help
```
```
Welcome to memod-s
usage: memod-s [-h] -i INPUT_DIRECTORY [-o OUTPUT_DIRECTORY] [-ml FILTLONG_MIN_LENGTH] [-kp FILTLONG_KEEP_PERCENT] [-rr RACON_ROUNDS]
                [-mr MEDAKA_ROUNDS] [-mm MEDAKA_MODEL] [-ed EGGNOG_DB] [-ab] [-ab_dir ABRICATE_DB_DIR] [-ab_db ABRICATE_DB_NAME]
                [-dm DORADO_MODELS [DORADO_MODELS ...]] [-dp DORADO_PU] [-q]

Snakefile wrapper for memod-s. For more details visit: https://github.com/AlessiaMarotta/memod-s

options:
  -h, --help            show this help message and exit
  -i INPUT_DIRECTORY, --input_directory INPUT_DIRECTORY
                        input directory containing FAST5 or POD5 files. (default: None)
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        output directory. (default: ./memod)

filtering options:
  -ml FILTLONG_MIN_LENGTH, --filtlong_min_length FILTLONG_MIN_LENGTH
                        minimum read length for filtlong. (default: 1000)
  -kp FILTLONG_KEEP_PERCENT, --filtlong_keep_percent FILTLONG_KEEP_PERCENT
                        percentage of best reads to keep. (default: 90)

assembly options:
  -rr RACON_ROUNDS, --racon_rounds RACON_ROUNDS
                        number of racon polishing rounds. (default: 4)
  -mr MEDAKA_ROUNDS, --medaka_rounds MEDAKA_ROUNDS
                        number of medaka polishing rounds. (default: 1)
  -mm MEDAKA_MODEL, --medaka_model MEDAKA_MODEL
                        medaka model to use. (default: r1041_e82_400bps_sup_v4.3.0)

annotation options:
  -ed EGGNOG_DB, --eggnog_db EGGNOG_DB
                        path to eggNOG database. (default: eggnog_db)
  -ab, --abricate       use abricate for resistance gene analysis. (default: False)
  -ab_dir ABRICATE_DB_DIR, --abricate_db_dir ABRICATE_DB_DIR
                        directory for abricate databases. (default: abricate_db_dir)
  -ab_db ABRICATE_DB_NAME, --abricate_db_name ABRICATE_DB_NAME
                        abricate database name. (default: vfdb)

dorado models:
  -dm DORADO_MODELS [DORADO_MODELS ...], --dorado_models DORADO_MODELS [DORADO_MODELS ...]
                        list of dorado models for basecalling. (default: ['dna_r10.4.1_e8.2_400bps_sup@v4.1.0',
                        'dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1', 'dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1'])
  -dp DORADO_PU, --dorado_pu DORADO_PU
                         If you don't have a GPU available, enter "cpu" as the processing unit to use for Dorado. (default: "cuda:all")

general options:
  -q, --quiet           suppress non-essential output. (default: False)


```
🛠 Basic usage, to run `memod-s` with default parameters:

```
./memod-s -i /path/to/your/input_dir/with/fast5/or/pod5
```
🔬 More advanced usage:

```
./memod-s --input_directory /path/to/your/input_dir/with/fast5/or/pod5 \
        --output_directory /path/to/output_dir \
        --filtlong_min_length 1500 \
        --filtlong_keep_percent 80 \
        --racon_rounds 5 \
        --medaka_rounds 2 \
        --medaka_model r1041_e82_400bps_sup_v4.3.0 \
        --dorado_models dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1 \
            dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1 \
        --ab ---ab_dir /path/to/abricate/db -ab_db vfdb
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
13. Restriction-Modification systems annotation with [MicrobeMod](https://github.com/cultivarium/MicrobeMod)
14. Combine genome-wide methylation profiles with genomic features with [MeStudio](https://github.com/combogenomics/MeStudio) 
15. Circular density plot for each motif with [Circlize](https://github.com/jokergoo/circlize)
16. Gene Set Enrichment Analysis with [fgsea](https://github.com/alserglab/fgsea)

### Bonus

1.  Antimicrobial resistance or virulence genes screening with [Abricate](https://github.com/tseemann/abricate)

## 💡 Results interpretation
### Methylation output
The two primary processed output files will be two tab-separated tables, one describing information for all methylated sites and one describing output for methylated motifs.
The output file `motifs.tsv` describe the consensus motif sequences called (Motif), the raw motif called by [STREME](https://meme-suite.org/meme/doc/streme.html) before cleaning up by MicrobeMod (Motif raw), the methylation type (Methylation type), the number of occurrence of the motif in the genome (Genome sites), the number of times the motif was methylated in the genome (Methylated sites), the ratio between genome sites and methylated sites (Methylation coverage), the mean of the percent of mapped reads that were methylated at all sites of the motif, with higher percentages indicating higher confidence in methylation calls (Average Percent Methylation per site). It might look something like:
| Index | Motif             | Motif_raw             | Methylation_type | Genome_sites | Methylated_sites | Methylation_coverage | Average_Percent_Methylation_per_site | Methylated_position_1 | Methylated_position_1_percent | Methylated_position_2 | Methylated_position_2_percent |
|-------|-------------------|-----------------------|------------------|--------------|------------------|----------------------|--------------------------------------|------------------------|-------------------------------|------------------------|-------------------------------|
| 0     | YGATCR            | 1-NNDNBYGATCRVNHNN    | 6mA              | 10690        | 9668             | 0.904                | 0.78                                 | 3                      | 50.0                          | 4                      | 50.0                          |
| 1     | No Motif Assigned | NA                    | 6mA              | NA           | 11373            | NA                   | 0.75                                 | NA                     | NA                            | NA                     | NA                            |
| 2     | CGSCG             | 1-NNHWNCGSCGNWDNN     | 5mC              | 8268         | 8216             | 0.994                | 0.82                                 | 2                      | 48.685                        | 4                      | 48.685                        |
| 3     | No Motif Assigned | NA                    | 5mC              | NA           | 279              | NA                   | 0.73                                 | NA                     | NA                            | NA                     | NA                            |

The `methylated_sites.tsv` output file is extensive and includes data on every methylated genomic site.

[MeStudio](https://github.com/combogenomics/MeStudio) generates a BED file for each feature, which includes (i) the seqid column indicating the name of each contig, (ii) the start and (iii) end positions of the feature, (iv) the gene ID present within the interval, (v) the number of methylations detected for that gene ID, and (vi) the corresponding protein product. 

Circular density plots are generated for the motifs across the genome, providing a visualization of methylation density in different genomic contexts. The outer circle of the plot displayed the genome annotation of the contigs, while each inner circle represents the different categories of methylated sites (CDS, nCDS, tIG and US). The plots can reveal distinct methylation patterns for each category within each motif. 

For example, in the circular plots below, upstream regions displayed a higher density of methylation compared to tIG regions, indicating a potential role of methylation in promoter regions:

<table>
  <tr>
    <td align="center">
      <img src="https://github.com/user-attachments/assets/3f88e676-9f01-484e-94ab-0f18393bc8db" width="300"/><br/>
      <sub><b>Motif A</b></sub>
    </td>
    <td align="center">
      <img src="https://github.com/user-attachments/assets/05b0c8cf-48c8-4cdd-ba21-08939094f1a9" width="300"/><br/>
      <sub><b>Motif B</b></sub>
    </td>
  </tr>
</table>


## 🗂️ Test memod-s

* You can use the sample data available at this [link](https://zenodo.org/records/15586873?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImM0NmU4YmVmLTQyYmUtNGQwZC1iOTAxLWI5MDA4NGYyYTc1ZCIsImRhdGEiOnt9LCJyYW5kb20iOiJiOGE0YTM3ZTg3NzU3ZjhhZTI1ZGEzMTZiM2Q3NDk4OSJ9.JOW0TFS_n7Imwjdel5Iffidzh6-a_FdT49jOQxdiUDuSvv0Sq9FcC1nTCSk69nAb1JTWwjpPl8Dup6vINTB8OQ) to test the workflow. This dataset contains a POD5 file from Nanopore sequencing of Vibrio aestuarianus. The file serves as an example input for users who wish to test the Snakemake-based memod-s workflow for methylation calling and analysis from Nanopore sequencing data.

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
