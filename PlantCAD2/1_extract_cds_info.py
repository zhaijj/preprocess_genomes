#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 13:56 2024/10/02
@author: JINGJING ZHAI (jz963@cornell.edu; zhaijingjing603@gmail.com)
Description: Extract CDS region coordinates, gene IDs, and transcript IDs from a GFF3 file.
"""

import pandas as pd
from collections import defaultdict
import sys
from tqdm import tqdm
import argparse

def extract_cds_coordinates(gff3_file):
    """
    Extract CDS coordinates (TSS to stop codon) and strand for each transcript from a GFF3 file.
    
    Parameters:
        gff3_file (str): Path to the GFF3 file.
    
    Returns:
        pd.DataFrame: DataFrame with columns ['Transcript ID', 'Chr', 'CDS Start', 'CDS End', 'Strand'].
    """
    transcript_data = []

    gff3_df = pd.read_csv(gff3_file, sep="\t", comment="#", header=None, 
                          names=["seqid", "source", "type", "start", "end", "score", 
                                 "strand", "phase", "attributes"])
    cds_df = gff3_df[gff3_df['type'] == 'CDS']
    
    cds_ranges = defaultdict(lambda: {"start": None, "end": None, "strand": None})

    # Extract CDS range and strand information
    for _, row in tqdm(cds_df.iterrows(), total=len(cds_df), desc="Extracting CDS info"):
        attributes = row['attributes']
        seqid = row['seqid']
        start = row['start']
        end = row['end']
        strand = row['strand']
        
        # Extract transcript ID (Parent attribute)
        transcript_id = None
        for attribute in attributes.split(";"):
            if attribute.startswith("Parent="):
                transcript_id = attribute.split("=")[1]
                break

        if transcript_id:
            # Update or initialize the CDS range for this transcript
            current_range = cds_ranges[transcript_id]
            current_range["seqid"] = seqid
            current_range["start"] = min(current_range["start"], start) if current_range["start"] is not None else start
            current_range["end"] = max(current_range["end"], end) if current_range["end"] is not None else end
            current_range["strand"] = strand

    # Convert the accumulated ranges to a list format
    for transcript_id, info in cds_ranges.items():
        transcript_data.append([transcript_id, info["seqid"], info["start"], info["end"], info["strand"]])
    
    return pd.DataFrame(transcript_data, columns=["Transcript ID", "Chr", "CDS Start", "CDS End", "Strand"])


def extract_gene_and_transcript_ids(gff3_file):
    """
    Extract transcript ID and gene ID from a GFF3 file.
    
    Parameters:
        gff3_file (str): Path to the GFF3 file.
    
    Returns:
        pd.DataFrame: DataFrame with columns ['Gene ID', 'Transcript ID'].
    """
    # List to hold rows for the DataFrame
    data = []

    # Read the GFF3 file into a pandas DataFrame
    gff3_df = pd.read_csv(gff3_file, sep="\t", comment="#", header=None, 
                          names=["seqid", "source", "type", "start", "end", "score", 
                                 "strand", "phase", "attributes"])
    
    # Filter rows where the type is mRNA (for transcript IDs)
    mRNA_df = gff3_df[gff3_df['type'] == 'mRNA']
    
    # Extract gene and transcript IDs
    for _, row in tqdm(mRNA_df.iterrows(), total=len(mRNA_df), desc="Extracting Gene/Transcript info"):
        attributes = row['attributes']
        
        transcript_id = None
        gene_id = None
        
        for attribute in attributes.split(";"):
            if attribute.startswith("ID="):
                transcript_id = attribute.split("=")[1]  # Extract transcript ID
            elif attribute.startswith("Parent="):
                gene_id = attribute.split("=")[1]  # Extract gene ID
        
        # Append to the data list if both transcript_id and gene_id are found
        if transcript_id and gene_id:
            data.append([gene_id, transcript_id])
    
    return pd.DataFrame(data, columns=["Gene ID", "Transcript ID"])


def main(gff3_file, output_file, mergeByGene=False, mode='max-CDS'):
    """
    Main function to extract and merge CDS coordinates and gene-transcript mappings.
    
    Parameters:
        gff3_file (str): Path to the GFF3 file.
        output_file (str): Path to save the output CSV file.
    """
    # Extract CDS coordinates and transcript-gene mappings
    cds_df = extract_cds_coordinates(gff3_file)
    gene_tx_df = extract_gene_and_transcript_ids(gff3_file)
    merged_df = pd.merge(gene_tx_df, cds_df, on="Transcript ID")

    if mergeByGene: 
        if mode == 'max-CDS':
            result = merged_df.groupby('Gene ID').agg({'Chr': 'first', 'CDS Start': 'min', 'CDS End': 'max', 'Strand': 'first'}).reset_index()
        elif mode == 'min-CDS':
            result = merged_df.groupby('Gene ID').agg({'Chr':'first', 'CDS Start': 'max', 'CDS End': 'min', 'Strand': 'first'}).reset_index()
        else:
            raise ValueError(f"Invalid merging mode: {mode}")
        result.to_csv(output_file, index=False, sep="\t")
        print(f"CDS information successfully saved to {output_file}")
    else:
        merged_df.to_csv(output_file, index=False, sep="\t")
        print(f"CDS information successfully saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract CDS region coordinates, gene IDs, and transcript IDs from a GFF3 file.")
    parser.add_argument("--gff3_file", type=str, help="Path to the GFF3 file.")
    parser.add_argument("--output_file", type=str, help="Path to save the output CSV file.")
    parser.add_argument("--mergeByGene", action='store_true', help="If set, the output will be merged by gene ID.")
    parser.add_argument("--mode", type=str, default='min-CDS', help="Merging mode for CDS coordinates. Options: 'max-CDS', 'min-CDS'. Default: 'min-CDS'.")
    args = parser.parse_args()
    
    gff3_file = args.gff3_file
    output_file = args.output_file
    mergeByGene = args.mergeByGene
    mode = args.mode
    
    main(gff3_file, output_file, mergeByGene, mode)