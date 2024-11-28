import logging
import pandas as pd
import numpy as np
import gzip
import sys

def get_sample_data(expression_data_path, samples, enst_list):
    sample_idx = None
    expression_rows = []  # List to collect rows

    with gzip.open(expression_data_path, 'rt') as gz_file:  # 'rt' mode for text
        for i, line in enumerate(gz_file):
            if i % 5000 == 0:
                logging.info(f"Processed {expression_data_path} to line {i}")
                logging.info(f"Size of loaded data: {sys.getsizeof(expression_rows)}")
            
            if sample_idx is None:
                # Parse the header (first line)
                pancan_samples = line.strip().split('\t')
                sample_idx = np.where(np.isin(pancan_samples, samples))[0].tolist()  # NumPy optimized index finding
                sample_idx.insert(0, 0)  # Add index for transcript identifier

                # Collect the column names (include transcript identifier)
                brca_pancan_samples = [pancan_samples[i] for i in sample_idx]
            else:
                # Process the data line by line
                pancan_count_data = line.strip().split('\t')
                enst_withversion = pancan_count_data[0]
                enst_withoutversion = enst_withversion.split('.')[0]
                if enst_withoutversion in enst_list:
                    pancan_count_data[0] = enst_withoutversion
                    pancan_count_data = np.array(pancan_count_data)  # Convert to NumPy array
                    sample_expression_data = pancan_count_data[sample_idx]  # Efficient indexing with NumPy
                    expression_rows.append(sample_expression_data)

    # Convert rows into a DataFrame once all rows are collected
    expression_df = pd.DataFrame(expression_rows, columns=brca_pancan_samples)
    
    return expression_df

def subset_common_columns(df1, df2):
    # Find the common columns between the two DataFrames
    common_columns = df1.columns.intersection(df2.columns)
    
    # Subset both DataFrames to only include the common columns
    df1_subset = df1[common_columns]
    df2_subset = df2[common_columns]
    
    # Return the altered DataFrames
    return df1_subset, df2_subset

def main():
    logging.basicConfig(level=logging.INFO)

    # Init paths
    counts_path = "data/raw/TcgaTargetGtex_expected_count.gz"
    tpm_path = "data/raw/TcgaTargetGtex_rsem_isoform_tpm.gz"
    counts_outpath = "data/processed/gtex_counts_subset.tsv"
    tpm_outpath = "data/processed/gtex_tpm_subset.tsv"
    deeploc_enst_path = "data/raw/cell_membrane_proteins_enst.txt"

    # Load in deeploc ENST
    deeploc_data = pd.read_csv(deeploc_enst_path, sep='\t')
    deeploc_ensts = set(deeploc_data["enst"].tolist())


    # Read in TCGA GTEX metadata
    meta_path = "data/raw/TCGA_GTEX_category.txt"
    metadata = pd.read_csv(meta_path, sep='\t')

    BRCA_meta = metadata[metadata["TCGA_GTEX_main_category"] == "GTEX Breast"]
    BRCA_samples = BRCA_meta["sample"].tolist()

    counts_df = get_sample_data(counts_path, BRCA_samples, deeploc_ensts)
    tpm_df = get_sample_data(tpm_path, BRCA_samples, deeploc_ensts)

    counts_df, tpm_df = subset_common_columns(counts_df, tpm_df)

    counts_df.to_csv(counts_outpath, sep="\t", index=False)
    logging.info(f"Counts subset written to {counts_outpath}")

    tpm_df.to_csv(tpm_outpath, sep="\t", index=False)
    logging.info(f"TPM subset written to {tpm_outpath}")

if __name__ == "__main__":
    main()