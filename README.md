# PyDESeq

## Introduction
The goal of this project is to build a tool in Python that calculates the log2 fold change values and p-values from RESM and TPM results to perform differential expression analysis of two datasets, control vs. treatment. Additionally, the tool will generate a volcano plot for differential expression with log2 fold change values as the x-axis and the -log10 p-value as the y-axis and prints out the top 10 most differentially expressed genes. 

[Prerequisites](#prerequisites) | [Install Instructions](#install) | [Basic Usage](#usage) | [pydeseq Options](#options) | [File formats](#format) | [Contributors](#credit)
<a name="prerequisites"></a>
## Prerequisites
`pydeseq` requires:
- Python 3.xx
- Python packages:
  - pandas
  - numpy
  - scipy
  - matplotlib

If these packages are not yet installed, use the `pip` command:
```
pip install pandas numpy scipy matplotlib
```
If you do not have root access, you can run the command above with the additional `--user` option to install locally:  
```
pip install --user pandas numpy scipy matplotlib
```

<a name="install"></a>
## Install Instructions

`pydeseq` can be installed with the following commands:
```
git clone https://github.com/shaniabuu/PyDESeq.git
cd PyDESeq
pip install .
```
Note: if you do not have root access, you can run the command above with additional options to install locally:
```
pip install . --user
```
To run `pydeseq` sucessfully, you working directory should be in the PyDESeq directory, otherwise, you might need to add this directory to your PATH by running:
```
export PATH=$PATH:~/.local/bin
```
Typing `pydeseq --help` should show a useful message and can be run to see if the install was successfull. 

<a name="usage"></a>
## Basic Usage
The basic usage of `./pydeseq.py` to process control and treatment files to calculate differential gene expression is:
```
pydeseq [-h] -c {CONTROLS ...} -t {TREATMENTS ...} -g {GENE_NAME_FILE} -p {PVALUE_THRESHOLD} [-o {OUTPUT_DIR}]
```
### pydeseq options
The required inputs for `pydeseq` are:
- `-c`, `--controls`: Control file. See [File Formats](#file-formats) for more details.
- `-t`, `--treatments`: Treatment file. See [File Formats](#file-formats) for more details.
- `-g`, `--gene_name_file`: Gene name file, for converting gene IDs to gene names. See [File Formats](#file-formats) for more details.
- `-p`, `--pvalue_threshold`: Sets p-value threshold for the output of the number of differentially expressed genes and the volcano plot.

Other additional options are:
- `-h`, `--help`: Shows the help message and exits
- `-o`, `--output_dir`: Sets output directy for output files. If not directory is given, saves output files to current working directory. 

<a name="format"></a>
## Input file formats

### Control and Treatment files
Control and Treatment files that are used should contain just two tab-delimited columns, the first for gene ids and the second for FPKM values:
```
gene_id      FPKM
```
column 1: `gene_id` stores the gene ids for the sequenced genes.

column 2: `FPKM` stores the FPKM (fragments per kilobase of transcript per million mapped reads) values of the sequenced genes.

Note: There should be two or more files for control and treatment files as replicates for the calculation of p-values.
See text files in [a link](https://github.com/shaniabuu/PyDESeq/tree/main/data/lab_data) for example.

### Gene name file
The gene name file that is used should contain just two tab-delimited columns, the first for gene ids and the second for gene names:
```
gene_id      gene_names
```
column 1: `gene_id` that includes those in control and treatment files.

column 2: `gene_names` for the readable gene names for each gene id.
See text files in [a link](https://github.com/shaniabuu/PyDESeq/tree/main/data/lab_data) for example.


Along with the files `differential_expression_results.csv` and `volcano_plot.png`. See [Output Files](#output-files) in [File Formats](#file-formats) for more details

<a name="output"></a>
## Output files
### `differential_expression_results.csv`
A comma-separated value (.csv) file with the differential expression analysis results stored in the following columns:
```
gene_id      log2FoldChange    pvalue
```
column 1: `gene_id` stores the gene ids for the identified differentially expressed genes.

column 2: `log2FoldChange` stores the values of the log2 fold change in expression of the gene between control(s) and treatment(s).

column 3: `p-value` stores the p-value for each gene where the null hypothesis is that that the gene expression is the same in control(s) vs. treatment(s). 

#### `volcano_plot.png`
A png image of the generated volcano plot with the inputted p-value threshold. A volcano plot is only generated if a p-value threshold is given. The red dots represent the differentially expressed genes with a p-value less than the inputted threshold, while the black dots represent the non-differentially expressed genes. The 10 genes are labeled with their gene names. 


<a name="credit"></a>
## Contributors
This repositiory was generated by Shania Bu, Amber Tse, and Janice Wu for our S24 CSE 185 final project, inspired by the R package [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). 
