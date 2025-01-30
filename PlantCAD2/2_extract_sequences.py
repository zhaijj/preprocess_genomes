#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10:48 2024/10/04
@author: JINGJING ZHAI (jz963@cornell.edu; zhaijingjing603@gmail.com)
Description: Extract psequences from the genome.
"""

import pandas as pd
import pysam
from tqdm import tqdm
import argparse
from Bio.Seq import Seq

def extract_sequences(genome_file, cds_info, flank_length = 1000):
    """
    Extract promoter and terminator sequences from the genome.
    
    Parameters:
        genome_file (str): Path to the genome file.
        cds_info (pd.DataFrame): DataFrame with columns ['Gene ID', 'Chr', 'CDS Start', 'CDS End', 'Strand'].
        flank_length (int): Length of the upstream and downstream regions to extract.
    
    Returns:
        pd.DataFrame: DataFrame with columns ['Gene ID', 'chrom', 'start', 'end', 'strand', 'Sequences'].
    """
    fasta = pysam.FastaFile(genome_file)
    cds_info = pd.read_csv(cds_info, sep="\t")

    results = []
    
    # Iterate over each transcript/gene and extract promoter/terminator sequences
    for _, row in tqdm(cds_info.iterrows(), total=len(cds_info), desc="Extracting flanking sequences"):
        gene_id = row['Gene ID']
        chrom = str(row['Chr'])
        cds_start = row['CDS Start'] - 1 - flank_length  # Convert to 0-based index
        cds_end = row['CDS End'] + flank_length
        strand = row['Strand']
        
        if chrom not in fasta.references:
            print(f"Chromosome {chrom} not found in the genome file.")
            continue

        # Get the length of the chromosome to avoid fetching sequences out of bounds
        chrom_length = fasta.get_reference_length(chrom)
        cds_start = max(0, cds_start)
        cds_end = min(chrom_length, cds_end)
        seq = fasta.fetch(reference=chrom, start=cds_start, end=cds_end)

        if strand == "-":
            seq = str(Seq(seq).reverse_complement())
        
        results.append({
            'Gene ID': gene_id,
            'chrom': chrom,
            'start': cds_start,
            'end': cds_end,
            'Strand': strand,
            'Sequences': seq,
        })
    
    # Close the fasta file to free resources
    fasta.close()
    
    return pd.DataFrame(results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract promoter and terminator sequences from the genome.")
    parser.add_argument("--genome_file", help="Path to the genome file.")
    parser.add_argument("--cds_info", help="Path to the CDS information file.")
    parser.add_argument("--flank_length", type=int, help="Length of the upstream region to extract.")
    parser.add_argument("--output_file", help="Path to save the output TSV file.")

    args = parser.parse_args()
    result_df = extract_sequences(args.genome_file, args.cds_info, args.flank_length)
    result_df.to_csv(args.output_file, sep="\t", index=False)