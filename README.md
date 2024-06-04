# PyDESeq

## Introduction
PyDESeq is a package to process control and treatment files to calculate differential gene expression. It calculates the log2 fold change and p-values for each gene and generates a filtered volcano plot. The script requires control and treatment FPKM (Fragments Per Kilobase Million) files and a gene name file.

[Prerequisites](#prerequisites) | [Install Instructions](#install) | [Basic Usage](#usage) | [Input file format](#format) | [Output files](#output) | [Example and testing](#example) | [Contributors](#credit)
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
pip install -e .
```
If you do not have root access, you can run the `pip install` command above with an additional option to install locally:
```
pip install --user -e .
```
To run `pydeseq` sucessfully, your working directory should be in the PyDESeq directory, otherwise, you might need to add this directory to your PATH. Note, this will have to be run again after opening a new terminal:
```
export PATH=$PATH:~/.local/bin
```
Typing `pydeseq --help` should show a useful message and can be run to see if the install was successfull. 

<a name="usage"></a>
## Basic Usage
The basic usage of `pydeseq` to process control and treatment files to calculate differential gene expression is:
```
pydeseq [-h] -c {CONTROLS ...} -t {TREATMENTS ...} -g {GENE_NAME_FILE} -p {PVALUE_THRESHOLD} [-o {OUTPUT_DIR}]
```
### pydeseq options
The required inputs for `pydeseq` are:
- `-c`, `--controls`: Control file. See [Input File Formats](#format) for more details.
- `-t`, `--treatments`: Treatment file. See [Input File Formats](#format) for more details.
- `-g`, `--gene_name_file`: Gene name file, for converting gene IDs to gene names. See [Input File Formats](#format) for more details.
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
- column 1: `gene_id` stores the gene ids for the sequenced genes.

- column 2: `FPKM` stores the FPKM (fragments per kilobase of transcript per million mapped reads) values of the sequenced genes.

Note: There should be two or more files for control and treatment files as replicates for the calculation of p-values.

Note: The file should contain a header for proper use. 

See text files in [/data/lab_data](https://github.com/shaniabuu/PyDESeq/tree/main/data/lab_data) for example.

### Gene name file
The gene name file that is used should contain just two tab-delimited columns, the first for gene ids and the second for gene names:

- column 1: `gene_id` that includes those in control and treatment files.

- column 2: `gene_names` for the readable gene names for each gene id.

Note: the file should not contain a header for proper use.

See GRCm38.75.gene_names in [/data/lab_data](https://github.com/shaniabuu/PyDESeq/tree/main/data/lab_data) for example.

<a name="output"></a>
## Output files
### text outputs
1. Top 10 Genes with the Smallest Non-zero P-values: The script prints the top 10 genes with the smallest non-zero p-values, showing their gene names, log2 fold changes, and p-values.
2. Number of Differentially Expressed Genes: The script prints the number of differentially expressed genes after filtering with the provided p-value threshold.

### output files
#### `differential_expression_results.csv`
A comma-separated value (.csv) file with the differential expression analysis results stored in the following columns:
```
gene_id      log2FoldChange    pvalue
```
column 1: `gene_id` stores the gene IDs for the identified differentially expressed genes.

column 2: `log2FoldChange` stores the values of the log2 fold change in expression of the gene between control(s) and treatment(s).

column 3: `p-value` stores the p-value for each gene where the null hypothesis is that that the gene expression is the same in control(s) vs. treatment(s). 

#### `volcano_plot.png`
This image file is a volcano plot of the log2 fold changes vs. -log10(p-values). Genes with p-values below the threshold are highlighted in red. The top 10 genes with the smallest non-zero p-values are labeled.

<a name="example"></a>
## Example and testing

- To test the package using [lab data](https://github.com/shaniabuu/PyDESeq/tree/main/data/lab_data):
```
cd PyDESeq
pydeseq -c data/lab_data/Chow_Rep1.txt data/lab_data/Chow_Rep2.txt  data/lab_data/Chow_Rep3.txt \
  -t data/lab_data/HFD_Rep1.txt data/lab_data/HFD_Rep2.txt data/lab_data/HFD_Rep3.txt \
  -g data/lab_data/GRCm38.75.gene_names  -p 0.05 -o test/lab_data/
```
The input files include three replices of mouse with standard "chow" diet as control and three replicates of mouse with high fat diet ("HFD") as treatment. The gene_id and gene_name conversion file is GRCm38.75.gene_names stored in the same directory as input files. The p-value threshold is set at 0.05 and the output file directory is test/lab_data/. 

See [/test/lab_data](https://github.com/shaniabuu/PyDESeq/tree/main/test/lab_data) for output with this command.

- To test the package using [public data](https://github.com/shaniabuu/PyDESeq/tree/main/data/public_data):
```
cd PyDESeq
pydeseq -c data/public_data/DOXY+_FUS-WT_REP1.txt data/public_data/DOXY+_FUS-WT_REP2.txt \
  -t data/public_data/DOXY+_FUS-P525L_REP1.txt data/public_data/DOXY+_FUS-P525L_REP2.txt \
  -g data/public_data/Human_gene_names.txt  -p 0.05 -o test/public_data/
```
```
cd PyDESeq
pydeseq -c data/public_data/DOXY-_FUS-WT_REP1.txt data/public_data/DOXY-_FUS-WT_REP2.txt \
  -t data/public_data/DOXY-_FUS-P525L_REP1.txt data/public_data/DOXY-_FUS-P525L_REP2.txt \
  -g data/public_data/Human_gene_names.txt  -p 0.05 -o test/public_data/
```

<a name="credit"></a>
## Contributors
This repositiory was generated by Shania Bu, Amber Tse, and Janice Wu for our S24 CSE 185 final project, inspired by the R package [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). 
