import pandas as pd
import argparse, os
from tqdm import tqdm
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

def split_sequence(row, window_size, step_size):
    """
    Split a sequence into windows of a given size and step size.
    """
    sequence = row['Sequences']
    sequence_splits = []
    
    # Split sequence into windows
    for start in range(0, len(sequence), step_size):
        window = sequence[start:(start + window_size)]
        if len(window) < window_size: 
            window = sequence[-window_size:]
            sequence_splits.append({**row, 'Sequences': window})
            break
            
        sequence_splits.append({**row, 'Sequences': window})
            
    return sequence_splits


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract promoter and terminator sequences from the genome.")
    parser.add_argument("--seq_dir", help="Path to the directory containing sequence files.")
    parser.add_argument("--window_size", type=int, default=512, help="Size of the window to split the sequences.")
    parser.add_argument("--step_size", type=int, default=256, help="Step size to move the window.")
    parser.add_argument("--threshold", type=float, default=0.1, help="Threshold for N content in sequences.")
    parser.add_argument("--output_file", help="Path to the output file.")
    parser.add_argument("--split", action='store_true', help="Whether to split train, validation and test sets.")

    args = parser.parse_args()

    seq_dir = args.seq_dir
    window_size = args.window_size
    step_size = args.step_size
    threshold = args.threshold
    # resDF = pd.DataFrame()
    # for file_name in os.listdir(seq_dir):
    #     if file_name.endswith('_genic_flank_5k.tsv'):
    #         logging.info(f"Processing file: {file_name}")
    #         seq_df = pd.read_csv(os.path.join(seq_dir, file_name), sep="\t")
    #         seq_df['assembly'] = file_name.split('_')[0]
    #         resDF = pd.concat([resDF, seq_df])
    
    # assembly_counts = resDF['assembly'].value_counts()
    # min_count = assembly_counts.min()
    # logging.info(f"Minimum count of sequences per assembly: {min_count}")
    
    # downsampled_df = resDF.groupby('assembly').apply(lambda x: x.sample(n=min_count, random_state=42)).reset_index(drop=True)
    # downsampled_df.to_csv('downsampled_gene.tsv', sep='\t', index=False)
    
    downsampled_df = pd.read_csv('downsampled_gene.tsv', sep='\t')
    logging.info(f"Processing {len(downsampled_df)} sequences")

    # remove sequences that are smaller than the window size
    downsampled_df = downsampled_df[downsampled_df['Sequences'].apply(lambda x: len(x) >= window_size)]

    df_expanded = downsampled_df.apply(lambda row: split_sequence(row, window_size, step_size), axis=1)

    logging.info("Expanding sequences")
    df_expanded = df_expanded.explode().reset_index(drop=True)
    df_expanded = df_expanded.dropna()
    df_expanded = pd.DataFrame(df_expanded.tolist())
    df_expanded = df_expanded[df_expanded['Sequences'].apply(lambda x: x.count('N') / len(x) <= threshold)]

    new_order = ['assembly', 'chrom', 'start', 'end', 'Strand', 'Sequences']
    df_expanded = df_expanded[new_order]
    new_names = {'Strand': 'strand', 'Sequences': 'seq'}
    df_expanded = df_expanded.rename(columns=new_names)


    df_expanded.to_csv(args.output_file, sep='\t', index=False)


    logging.info(f"Downsampled sequences saved to: {args.output_file}")

    if args.split:
        train_df = df_expanded.sample(frac=0.8, random_state=42)
        val_test_df = df_expanded.drop(train_df.index)
        val_df = val_test_df.sample(frac=0.5, random_state=42)
        test_df = val_test_df.drop(val_df.index)

        train_df.to_csv(args.output_file.replace('.tsv', '_train.tsv'), sep='\t', index=False)
        val_df.to_csv(args.output_file.replace('.tsv', '_val.tsv'), sep='\t', index=False)
        test_df.to_csv(args.output_file.replace('.tsv', '_test.tsv'), sep='\t', index=False)
        logging.info(f"Train, validation and test sets saved to: {args.output_file.replace('.tsv', '_train.tsv')}, {args.output_file.replace('.tsv', '_val.tsv')}, {args.output_file.replace('.tsv', '_test.tsv')}")
