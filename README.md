# Pre-process genomes before pre-training

## 1. chunk genomes into windows

### Tool requirements
- bedtools

### Dataset requirements
- softmasked genomes in fasta format
- genome index files in fasta.fai format

```bash
window_size=512
step_size=256

# create the output folder if it does not exist
if [ ! -d 1_windows ]; then
    mkdir 1_windows
fi

ls 1_softmasked_genomes/*fai | while read line
do
    echo ${line}
    prefix=`basename ${line} | sed 's/.fasta.fai//g'`
    bedtools makewindows -g ${line} -w ${window_size} -s ${step_size} > 1_windows/${prefix}.windows_512.bed
done
```

## 2. Assign each genomic window to a unique class

### R packages requirements
- rtracklayer
- snowfall
- GenomicFeatures
- Biostrings
- argparse
- parallel


### Datasets requirements
- Genome annotation in .GFF3 format
- RepeatMasker output files OR TE in .GFF3 format
- BED files of genomic windows

```bash
# 32 the number of cores to use
Rscript genome_window_class.R input.bed annotation.gff3 repeatmasker.gff3 output.bed 32 output.txt
```

## 3. Down sample the windows

Use the script `downSampleDataset.R` to down sample the windows, it's not open-boxed yet, but you can modify it to fit your needs.