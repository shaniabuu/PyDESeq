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

`pydeseq` requires the following python libraries to be installed:

pandas
numpy 
scipy
matplotlib


pydeseq can be installed with the following commands:
```
git clone https://github.com/shaniabuu/PyDESeq.git
cd PyDESeq/pydeseq

```

Typing `./pydeseq.py --help` should show a useful message and can be run to see if the install was successfull. 

<a name="usage"></a>
## Basic Usage
The basic usage of `./pydeseq.py` to process control and treatment files to calculate differential gene expression is:
```
./pydeseq.py [-h] -c {CONTROLS ...} -t {TREATMENTS ...} -g {m,h} [-p {PVALUE_THRESHOLD}] [-o {OUTPUT_DIR}]
```
To run `./pydeseq.py` on our test example files:
```
./pydeseq.py -c ../data/lab_data/Chow_Rep1.txt ../data/lab_data/Chow_Rep2.txt ../data/lab_data/Chow_Rep3.txt -t ../data/lab_data/HFD_Rep1.txt ../data/lab_data/HFD_Rep2.txt ../data/lab_data/HFD_Rep3.txt -g m -p 0.05
```

This should produce the following output:

```
                  gene_name  log2FoldChange        pvalue
gene_id
ENSMUSG00000081664   Gm15544      -15.609669  9.629650e-33
ENSMUSG00000048251    Bcl11b        0.807332  1.003089e-32
ENSMUSG00000034987      Hrh2       -0.678056  1.504633e-32
ENSMUSG00000017417    Plxdc1        0.321927  3.530872e-32
ENSMUSG00000026730      Pter        0.654322  1.116894e-04
ENSMUSG00000033174      Mgll        0.990863  1.729810e-04
ENSMUSG00000021136     Smoc1        0.152529  1.964658e-04
ENSMUSG00000020261   Slc36a1        0.423822  2.933470e-04
ENSMUSG00000038422     Hdhd3        0.728353  3.483287e-04
ENSMUSG00000032667      Pon2        0.220157  3.778721e-04
Number of differentially expressed genes (p-value > 0.05): 19692
```
Along with the files `differential_expression_results.csv` and `volcano_plot.png`. See [Output Files](#output-files) in [File Formats](#file-formats) for more details

<a name="options"></a>
## pydeseq options
The required inputs for `./pydeseq.py` are:
- `-c`, `--controls`: Control file(s). 
- `t`, `--treatments`: Treatment file(s).
- `-g`, `--genome`: Genome type. Use `m` for mouse, `h` for human

Other additional options are:
- `-h`, `--help`: Shows the help message and exits
- `p`, `--pvalue_threshold`: Sets p-value threshold for the volcano plot. If no p-value threshold is inputed, a volcano plot will not be generated. 
- `-o`, `--output_dir`: Sets output directy for output files. If not directory is given, saves output files to current working directory. 

<a name="format"></a>
## File formats

### Inupt files (from test data)
#### `Chow_Rep1.txt Chow_Rep2.txt Chow_Rep3.txt`

Text files for the replicates of the control samples of mice fed with a Chow diet. They contain two columns :
```
gene_id      FPKM
```
column 1: `gene_id` stores the gene ids for the sequenced genes

column 2: `FPKM` stores the FPKM (fragments per kilobase of transcript per million mapped reads) values of the sequenced genes.


#### `HFD_Rep1.txt HFD_Rep2.txt HFD_Rep3.txt`

Text files for the replicates of the treatment samples of mice fed with a high fat diet. They contain two columns :
```
gene_id      FPKM
```
column 1: `gene_id` stores the gene ids for the sequenced genes

column 2: `FPKM` stores the FPKM (fragments per kilobase of transcript per million mapped reads) values of the sequenced genes.

Control and Treatment files that are used should contain just two tab-delimited columns, the first for gene ids and the second for FPKM values. 

### Output files
#### `differential_expression_results.csv`
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
