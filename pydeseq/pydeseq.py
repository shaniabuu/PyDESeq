#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_rel
import os
import matplotlib.pyplot as plt

def load_gene_names(file_path):
    return pd.read_csv(file_path, sep='\t', header=None, names=['gene_id', 'gene_name'], index_col='gene_id')

def read_fpkm_file(file_path):
    return pd.read_csv(file_path, sep='\t')

def calculate_log2_fold_change(ctrl_avg, treat_avg, pseudocount=1e-6):
    return np.log2((treat_avg + pseudocount) / (ctrl_avg + pseudocount))

def main():
    parser = argparse.ArgumentParser(
        description="Process control and treatment files to calculate differential gene expression.",
        usage="pydeseq.py [-h] -c {CONTROLS ...} -t {TREATMENTS ...} -g {m,h} [-p {PVALUE_THRESHOLD}] [-o {OUTPUT_DIR}]"
    )
    parser.add_argument('-c', '--controls', nargs='+', required=True, help='Control file(s)')
    parser.add_argument('-t', '--treatments', nargs='+', required=True, help='Treatment file(s)')
    parser.add_argument('-g', '--genome', required=True, choices=['m', 'h'], help='Genome type: m for mouse, h for human')
    parser.add_argument('-p', '--pvalue_threshold', type=float, help='P-value threshold for volcano plot')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')

    args = parser.parse_args()

    # Load gene names
    gene_name_file = 'GRCm38.75.gene_names' if args.genome == 'm' else 'hg38.gene_names'
    gene_names = load_gene_names(os.path.join(os.path.dirname(__file__), gene_name_file))

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

    combined_df = combined_df.join(gene_names, on='gene_id')
    combined_df_nonzero_pvalue = combined_df[combined_df['pvalue'] > 0]

    # Print top 10 genes with the smallest non-zero p-values
    top_genes_nonzero_pvalue = combined_df_nonzero_pvalue.nsmallest(10, 'pvalue')
    print(top_genes_nonzero_pvalue[['gene_name', 'log2FoldChange', 'pvalue']])

    # Print number of differentially expressed genes
    num_diff_expr_genes = combined_df[combined_df['pvalue'] > 0.05].shape[0]
    print("Number of differentially expressed genes (p-value > 0.05):", num_diff_expr_genes)

    # Create volcano plot if p-value threshold is specified and p-value column exists
    if args.pvalue_threshold and 'pvalue' in combined_df.columns:
        # Add a pseudocount to p-values to avoid log10(0)
        pseudocount_pvalue = 1e-10
        combined_df['-log10(pvalue)'] = -np.log10(combined_df['pvalue'] + pseudocount_pvalue)
        plt.scatter(combined_df['log2FoldChange'], combined_df['-log10(pvalue)'], color='black')

        # Highlight significant points
        significant = combined_df['pvalue'] < args.pvalue_threshold
        plt.scatter(combined_df.loc[significant, 'log2FoldChange'],
                    combined_df.loc[significant, '-log10(pvalue)'],
                    color='red')

        # Add horizontal threshold line
        plt.axhline(-np.log10(args.pvalue_threshold), linestyle='--', color='gray')

        # Label top 10 genes
        if '-log10(pvalue)' in combined_df.columns:
            top_genes_to_label = top_genes_nonzero_pvalue.copy()
            top_genes_to_label['-log10(pvalue)'] = combined_df.loc[top_genes_nonzero_pvalue.index, '-log10(pvalue)']
            for _, row in top_genes_to_label.iterrows():
                plt.text(row['log2FoldChange'], row['-log10(pvalue)'], row['gene_name'])

        plt.xlabel('log2 Fold Change')
        plt.ylabel('-log10(p-value)')
        plt.title('Volcano Plot')
        # Adjust figure size
        plt.gcf().set_size_inches(12, 6)  # Set width to 12 inches and height to 6 inches
        plt.savefig(os.path.join(args.output_dir, 'volcano_plot.png'), bbox_inches='tight')

        
if __name__ == '__main__':
    main()
