# `memod-s`

> **Note**

> An easy-to-use workflow for generating BLAH BLAH

![Snakemake](https://img.shields.io/badge/Snakemake-%3E%3D5.10.0%2C%3C5.31.1-green?style=for-the-badge) ![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) ![C](https://img.shields.io/badge/c-%2300599C.svg?style=for-the-badge&logo=c&logoColor=white) ![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white) ![Bash](https://img.shields.io/badge/bash-%234EAA25.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

Hi :blush:

`memod-s` is a Snakemake workflow that BLAH BLAH

## ⚙️ Installation

You can start using `memod-s` on your cluster with just one line of code with the [mamba package manager](https://github.com/mamba-org/mamba)

```
mamba create -n memod -c bioconda memod
```

This will create an environment called `memod` and start installing dependencies.

## 🔧 Usage

Clone this repo

```
git clone https://github.com/alenana99/memod_snake.git
```

Run `memod-s` without any arguments to see usage instructions:

```
python memod-s.py
```
```
Usage: bash memod-s.py 
                       
 Options:
```

## 🐍 Workflow

### Core

1. Convert fast5 files to pod5 with [Pod5 File Format Documentation](https://pod5-file-format.readthedocs.io/en/latest/docs/tools.html#pod5-convert-fast5)
2. Download specific methylation models with Dorado. Check the [avaiable basecalling models](https://github.com/nanoporetech/dorado?tab=readme-ov-file#available-basecalling-models)
3. Basecalling with [Dorado](https://github.com/nanoporetech/dorado)
4. Map reads to reference with minimap2
5. Quality check with
6. Quality filtering with
7. Quality check after fitlering 
8. Assembly with dragonflye
9. Assembly evaluation with quast
10. Genome annotation with prokka
11. Functional annotation with eggnog
12. Call methylation with [MicrobeMod](https://github.com/cultivarium/MicrobeMod)
13. MeStudio 
14. Circular density plot for each motif with 
15. Gene Set Enrichment Analysis with fgsea

### Bonus

1.  Abricate 

## 🖇️ Publications

The `memod-s` workflow BLAH BLAH

```
NOSTRO ARTICOLO
```
```
