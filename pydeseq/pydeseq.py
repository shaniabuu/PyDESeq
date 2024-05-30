#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_rel
import os
import matplotlib.pyplot as plt
import warnings

# Suppress specific RuntimeWarning from scipy.stats
warnings.filterwarnings("ignore", category=RuntimeWarning, message="Precision loss occurred in moment calculation due to catastrophic cancellation. This occurs when the data are nearly identical. Results may be unreliable.")

def load_gene_names(file_path):
    return pd.read_csv(file_path, sep='\t', header=None, names=['gene_id', 'gene_name'], index_col='gene_id')

def read_fpkm_file(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t')
        if 'gene_id' not in df.columns or 'FPKM' not in df.columns:
            raise ValueError(f"File {file_path} does not contain required columns 'gene_id' and 'FPKM'.")
        return df
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        exit(1)

def calculate_log2_fold_change(ctrl_avg, treat_avg, pseudocount=1e-6):
    return np.log2((treat_avg + pseudocount) / (ctrl_avg + pseudocount))

class MyArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f'Error: {message}\n')
        self.print_help()
        exit(2)

def main():
    parser = MyArgumentParser(
        description="Process control and treatment files to calculate differential gene expression.",
        usage="pydeseq [-h] -c {CONTROLS ...} -t {TREATMENTS ...} -g {GENE_NAME_FILE} -p {PVALUE_THRESHOLD} [-o {OUTPUT_DIR}]"
    )
    parser.add_argument('-c', '--controls', nargs='+', required=True, help='Control files')
    parser.add_argument('-t', '--treatments', nargs='+', required=True, help='Treatment files')
    parser.add_argument('-g', '--gene_name_file', required=True, help='Gene name file')
    parser.add_argument('-p', '--pvalue_threshold', type=float, required=True, help='P-value threshold for defining significance')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')

    args = parser.parse_args()

    # Load gene names
    gene_names = load_gene_names(args.gene_name_file)

    # Load and average FPKM values
    control_data = [read_fpkm_file(file) for file in args.controls]
    treatment_data = [read_fpkm_file(file) for file in args.treatments]

    control_df = pd.concat(control_data).groupby('gene_id').mean()
    treatment_df = pd.concat(treatment_data).groupby('gene_id').mean()

    # Calculate log2 fold change with pseudocount
    combined_df = pd.DataFrame({
        'ctrl_avg': control_df['FPKM'],
        'treat_avg': treatment_df['FPKM']
    })
    combined_df['log2FoldChange'] = calculate_log2_fold_change(combined_df['ctrl_avg'], combined_df['treat_avg'])

    # Perform paired t-test
    all_ctrl_fpkm = pd.concat([df.set_index('gene_id')['FPKM'] for df in control_data], axis=1)
    all_treat_fpkm = pd.concat([df.set_index('gene_id')['FPKM'] for df in treatment_data], axis=1)

    t_stat, p_values = ttest_rel(all_ctrl_fpkm, all_treat_fpkm, axis=1)
    combined_df['pvalue'] = p_values

    # Save output file
    os.makedirs(args.output_dir, exist_ok=True)
    output_file = os.path.join(args.output_dir, 'differential_expression_results.csv')
    combined_df.to_csv(output_file, columns=['log2FoldChange', 'pvalue'])
    
    # Join gene names
    combined_df = combined_df.join(gene_names, on='gene_id')
    combined_df_nonzero_pvalue = combined_df[combined_df['pvalue'] > 0].copy()
    
    pseudocount_pvalue = 1e-10
    # Create -log10(pvalue) column
    combined_df_nonzero_pvalue.loc[:, '-log10(pvalue)'] = -np.log10(combined_df_nonzero_pvalue['pvalue'] + pseudocount_pvalue)

    # Filter out extreme log2 fold change and p-values
    filtered_combined_df_nonzero_pvalue = combined_df_nonzero_pvalue[
        (combined_df_nonzero_pvalue['log2FoldChange'].between(-10, 10)) &
        (combined_df_nonzero_pvalue['-log10(pvalue)'] <= 9)
    ]

    # Print filtered top 10 genes with the smallest non-zero p-values
    filtered_top_genes_nonzero_pvalue = filtered_combined_df_nonzero_pvalue.nsmallest(10, 'pvalue').reset_index(drop=True)
    print(filtered_top_genes_nonzero_pvalue[['gene_name', 'log2FoldChange', 'pvalue']])

    # Print number of differentially expressed genes after filtering
    num_filtered_diff_expr_genes = filtered_combined_df_nonzero_pvalue[filtered_combined_df_nonzero_pvalue['pvalue'] > args.pvalue_threshold].shape[0]
    print(f"Number of differentially expressed genes after filtering (p-value > {args.pvalue_threshold}):", num_filtered_diff_expr_genes)

    # Create filtered volcano plot if p-value threshold is specified and p-value column exists
    plt.scatter(filtered_combined_df_nonzero_pvalue['log2FoldChange'], filtered_combined_df_nonzero_pvalue['-log10(pvalue)'], color='black')

    # Highlight significant points
    significant = filtered_combined_df_nonzero_pvalue['pvalue'] < args.pvalue_threshold
    plt.scatter(filtered_combined_df_nonzero_pvalue.loc[significant, 'log2FoldChange'],
                filtered_combined_df_nonzero_pvalue.loc[significant, '-log10(pvalue)'],
                color='red')

    # Add horizontal threshold line
    plt.axhline(-np.log10(args.pvalue_threshold), linestyle='--', color='gray')

    # Label top 10 genes
    for _, row in filtered_top_genes_nonzero_pvalue.iterrows():
        plt.text(row['log2FoldChange'], row['-log10(pvalue)'], row['gene_name'])

    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10(p-value)')
    plt.title('Filtered Volcano Plot')
    # Adjust figure size
    plt.gcf().set_size_inches(12, 6)  # Set width to 12 inches and height to 6 inches
    plt.savefig(os.path.join(args.output_dir, 'volcano_plot.png'), bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main()
