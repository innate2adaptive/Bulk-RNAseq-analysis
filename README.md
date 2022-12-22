# Bulk-RNAseq-analysis
Repo with scripts used for bulk RNAseq analysis in the Noursadeghi's lab 

The repo contains:

- Zscore_gene_expression_module_analysis: all scripts necessary to run a Z-score module analysis with bulk RNAseq expression data. The folder contains the target gene modules (full_ipa.csv, created with Ingenuity Pathway Analysis) and 3 scripts for the analysis. The input file is a data frame with log2 transcript per million (TPM) per gene per sample.  The first step calculates the average gene-pair correlation in the gene modules and in 100 iterations of randomly selected gene clusters of the same sizes (CORRELATION_SCRIPT.R). Another script runs a similar analysis using average gene expression rather than correlation values (Avg_Expression.R). To enable comparison across data sets, a third script transforms average gene expression values for each module to standardized (Z scores) using mean and standard deviation of ram=randomly selected gene sets. Statistical significant differences in Z scores between groups is then identified by t-tests (with multiple testing correction). 


