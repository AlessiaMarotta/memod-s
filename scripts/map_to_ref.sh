#!/bin/bash

BAM=$1
REF=$2
OUT=$3

samtools fastq "$BAM" -T MM,ML | \
minimap2 -t 14 --secondary=no -ax map-ont -y "$REF" - | \
samtools view -b - | \
samtools sort -@10 -o "${OUT}" -

samtools index "${OUT}"
