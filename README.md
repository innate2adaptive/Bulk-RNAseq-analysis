# Bulk-RNAseq-analysis
Repo with scripts used for bulk RNAseq analysis in the Noursadeghi's lab 

The repo contains scripts to run cytokine analysis with bulk RNAseq data.
- Cytokine regulators of genes enriched in the tuberculine skin test (E-MTAB-6816) were identified using Ingenuity Pathway Analysis (IPA_GenesModules.txt)
- Correlation.R: this script takes a file with log2 transcript per million (TPM) per gene per sample and calculates the average gene-pair correlation. It then compares the average correlation of the cytokine modules (IPA_GenesModules.txt) and compares it with 100 iterations of randomly selected gene modules of the same sizes.
- Average.R: similar to Correlation.R, this script uses the average log2 TPM expression rather than correlation values
