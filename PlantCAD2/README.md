# PlantCAD2 Dataset Preprocessing

This document describes the preprocessing pipeline for the PlantCAD2 pre-training dataset.

## Overview

The pipeline consists of three main steps:
1. Genome download and repeat masking
2. CDS information extraction
3. Flanking sequence extraction

## 1. Genome Download and Repeat Masking

### 1.1 Data Collection

We collected 65 plant genomes from [Phytozome](https://phytozome-next.jgi.doe.gov/). To ensure taxonomic balance, we manually selected one representative genome per genus.

### 1.2 Data Organization

The dataset is organized in the following structure:

```
├── softmasked_genomes/
│   ├── Acoerulea.softmasked.fa
│   ├── Acomosus.softmasked.fa
│   └── Ahypochondriacus.softmasked.fa
│
├── annotation/
│   ├── Acoerulea.gff3
│   ├── Acomosus.gff3
│   └── Ahypochondriacus.gff3
```

### 1.3 De-novo Repeat Masking Pipeline

The following script performs genome-wide k-mer frequency-based repeat masking, we used a [pipeline](https://github.com/baoxingsong/dCNS) here: 

```bash
#!/bin/bash

function mask_genome() {
    input_fa="assemblies/$1"
    prefix=$(basename "${input_fa}" .fa)
    output_dir="softmasked_genomes"
    intermediate_dir="kmer_masking_intermediate_output"

    # Create output directories
    mkdir -p "${output_dir}" "${intermediate_dir}"

    # Generate k-mer histogram
    kat hist -t 4 "${input_fa}" -m 20 -o "${intermediate_dir}/${prefix}.kat.m20.hist"

    # Calculate frequency threshold
    Rscript getFreqThreshold.R "${intermediate_dir}/${prefix}.kat.m20.hist" \
        "${intermediate_dir}/${prefix}.kfreq.cutoff"

    # Count k-mers
    kat_jellyfish count -m 20 -s 100M -t 4 -C "${input_fa}" \
        -o "${intermediate_dir}/${prefix}_k20_count.js"

    # Export k-mers
    kat_jellyfish dump "${intermediate_dir}/${prefix}_k20_count.js" \
        > "${intermediate_dir}/${prefix}_k20_count_dumps.fa"

    # Perform genome masking
    /home/jz963/pixi_envs/biosoft/dCNS/dCNS maskGenome -d \
        -i "${input_fa}" \
        -o "${output_dir}/${prefix}.softmasked.fa" \
        -k "${intermediate_dir}/${prefix}_k20_count_dumps.fa" \
        -f $(cat "${intermediate_dir}/${prefix}.kfreq.cutoff")
}

export -f mask_genome

# Process genomes in parallel
cat genome_prefix.txt | parallel -j 8 mask_genome {}
```

## 2. CDS Information Extraction

Extract regions upstream of translation initiation sites and downstream of translation stop sites:

```bash
for i in annotation/*gff3; do
    prefix=$(basename $i .gff3)
    python 1_extract_cds_info.py \
        --gff3_file ${i} \
        --output_file ${prefix}_cds_info.tsv \
        --mergeByGene \
        --mode 'max-CDS'
done
```

## 3. Flanking Sequence Extraction

Extract 5kb flanking sequences for each gene:

```bash
ls softmasked_genomes/*fa | while read line; do
    prefix=$(basename ${line} .softmasked.fa)
    cds_info="annotation/${prefix}_cds_info.tsv"
    output_file="annotation/${prefix}_genic_flank_5k.tsv"
    
    python 2_extract_sequences.py \
        --genome_file "softmasked_genomes/${prefix}.softmasked.fa" \
        --cds_info ${cds_info} \
        --flank_length 5000 \
        --output_file ${output_file} &
done
```
