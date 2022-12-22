setwd("~/Documents/R/correlation_analysis_TST/data")

########## 1. Average cluster target gene co-correlation ##########
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)

#abortive infection data set - normalised counts

## A. Input files - replace as necessary

transcriptome<- read.csv("DeduplicatedLog2TPM_AbortiveInfection_FUP4_Deseq2FDR05_Up_2022-11-7.csv", header=T)
rownames(transcriptome) <- transcriptome[,1]
transcriptome <- transcriptome[,-1]
ipa <- read.delim("uraZ2_AbortiveInfectionFUP4vNegNeg_Up.txt",header=T) # IPA output file

## B. Create a co-correlation matrix of the transcriptome
transcriptome <- transcriptome[order(rownames(transcriptome)),]
corr_mat <- cor(t(transcriptome), method = "spearman")
corr_mat
corr_mat[is.na(corr_mat)] <- 0


## C. Format the IPA output file
types <- c("cytokine", "transcription regulator", "kinase", "transmembrane receptor")
colnames(ipa)
ipa <- ipa %>%
  filter(B.H.corrected.p.value < 0.05  & Molecule.Type %in% types) %>%
  select(Upstream.Regulator, Molecule.Type, B.H.corrected.p.value, Target.Molecules.in.Dataset)
colnames(ipa) <- c("regulator","type","FDR","target")
ipa <- tibble::rowid_to_column(ipa, "ID") #add an ID column corresponding to each cluster
ipa$target <- as.character(ipa$target) #target column is integer -> convert to character for next step
ipa_sep <- separate_rows(ipa, target, sep = ",|/", convert = TRUE) #separate clusters into individual rows
ipa_sep_reduced <- subset(ipa_sep, target %in% colnames(corr_mat)) #only the genes found in the transcriptome co-correlation matrix
ipa_sep_reduced %>%
  group_by(ID) %>%
  mutate("count" = n()) -> ipa_sep_reduced # add a column with the cluster size

## D. Calculate the average of all pairwise correlations for the target genes of each upstream regulator cluster
avg_interactomes <- matrix(nrow=0, ncol=2) # empty data frame to add the calculated average correlation values
colnames(avg_interactomes) <- c("regulator", "avg_corr")
for (i in 1:max(ipa_sep_reduced$ID)){
  ipa_sep_reduced %>%
    filter(ID == i) -> hub
  genes <- hub$target[which(hub$target %in% row.names(corr_mat))]
  corr_mat_network <- corr_mat[genes,genes]
  lower_matrix <- lower.tri(corr_mat_network)
  corr_mat_network_lower <- corr_mat_network[lower_matrix]
  mean_interactome <- mean(corr_mat_network_lower[corr_mat_network_lower>=0])
  avg_interactomes <- rbind(avg_interactomes, c(as.character(unique(hub$regulator)), mean_interactome))
}
# add the calculated cluster average correlations to a data frame containing the other information
ipa_sep_reduced$target <- NULL
ipa_sep_reduced <- unique(ipa_sep_reduced)
ipa_sep_reduced <- merge(ipa_sep_reduced, avg_interactomes, by.x = "regulator", by.y = "regulator") # merge the ipa output file and the correlation values by common regulator name


## E. Output file
setwd("~/Documents/R/correlation_analysis_TST/output")
write.csv(ipa_sep_reduced, "1_avg-correl_abortinf_amendTPM.csv", row.names = F)


ipa_sep_reduced <- read.csv("1_avg-correl_abortinf_nomcounts.csv")


# generate random distributions
#check max count
print(max(ipa_sep_reduced$count))
print(ipa_sep_reduced$count)
count <- ipa_sep_reduced$count

p <- 210
p <- count
r <- 1000
gene_names <- rownames(corr_mat)
gene_names <- rownames(transcriptome)
gene_names <- unique(gene_names)
gene_names
length(gene_names)

random_distributions_1000 <- matrix(nrow = 0, ncol = 6) # empty matrix to fill with the distribution values
colnames(random_distributions_1000) <- c("size", "mean", "sd", "84.13%", "97.72%",	"99.87%")
for (k in p){
  random_genes <- matrix(nrow = 0, ncol = 2)
  colnames(random_genes) <- c("iteration", "average_correlation")
  for (i in (1:r)){ # Do this 100 times (could also try doing 1000 times if it doesn't take too long for the loop to run)
    genes <- sample(gene_names, k) # random sample from transcriptome
    corr_mat_network <- corr_mat[genes,genes]
    lower_matrix <- lower.tri(corr_mat_network)
    corr_mat_network_lower <- corr_mat_network[lower_matrix]
    mean_interactome <- mean(corr_mat_network_lower[corr_mat_network_lower>=0]) # Calculate the average correlation for each of the 100 (or 1000) random clusters
    random_genes <- rbind(random_genes, c(i, mean_interactome))
  }
  mean <- mean(as.numeric(random_genes[,2]))
  sd <- sd(as.numeric(random_genes[,2]))
  random_distributions_1000 <- rbind(random_distributions_1000, c(k, mean, sd, mean + sd, mean + 2*sd, mean + 3*sd))
}
colnames(random_distributions_1000) <- c("size", "mean", "sd", "X84.13", "X97.72",	"X99.86")

## D. Output file
#write.csv(random_distributions_100,

write.csv(random_distributions_1000,
          paste0("random_distributions_",203,"size_",r,"abortiveinf_TPMamend",Sys.Date(),".csv"),
          row.names = F)

#check if you are not using ipa output to calculate random correlations - passed!


ipa_avg_corr <-read.csv("1_avg-correl_abortinf_amendTPM.csv") # Output file from Script 1

random_distributions_1000 <- read.csv("random_distributions_203size_1000abortiveinf_TPMamend2022-11-08.csv") # output file from Script 2 

## B. Determining which clusters are FDR significant
#ipa_avg_corr <- ipa_avg_corr[ipa_avg_corr$count >= 4,] # take out clusters smaller than 4
ipa_avg_corr$rc_mean <- random_distributions_1000$mean[match(ipa_avg_corr$count, random_distributions_1000$size)] #mean
ipa_avg_corr$rc_sd <- random_distributions_1000$sd[match(ipa_avg_corr$count, random_distributions_1000$size)] #sd
ipa_avg_corr$zscore <- (ipa_avg_corr$avg_corr - ipa_avg_corr$rc_mean)/ipa_avg_corr$rc_sd #zscores
ipa_avg_corr$pvalue <- pnorm(ipa_avg_corr$zscore, lower.tail = FALSE) #p-values
#ipa_less_0.05 <- ipa_avg_corr[which(ipa_avg_corr$pvalue <= 0.05),] # p-values <= 0.05
fdr_threshold <- nrow(ipa_avg_corr)
ipa_avg_corr$adj_pvalue <- ipa_avg_corr$pvalue*fdr_threshold # calculating the adjusted p-value
ipa_significant <- ipa_avg_corr[which(ipa_avg_corr$adj_pvalue <= 0.05),] # selecting clusters with adj. p-value <= 0.05

## C. Output files
#write.csv(ipa_significant, "../data/UR_module_identification/3_fdrsig_40size_inf_pre0.csv", row.names = F) # clusters with adj. p-value <= 0.05
#write.csv(ipa_significant, "script3_significantcorrelations_abortinf_nomcounts.csv", row.names = F) # clusters with adj. p-value <= 0.05
write.csv(ipa_avg_corr, "script3_ipaavgcorr_abortinf_TPMamend.csv", row.names = F)
## D. Plotting the FDR-significant clusters
library(ggplot2)

# plotting all the clusters against the frequency distribution and making FDR-sig clusters blue
ipa_avg_corr$Significance <- (ipa_avg_corr$regulator%in%ipa_significant$regulator) # need to have a T/Fsignificance column
colnames(random_distributions_1000)

random_genes_plot_coloursig <- ggplot(ipa_avg_corr, aes(avg_corr, count, colour = Significance)) + 
  xlab("Average correlation coefficient") + ylab("Network size (number of target genes)") + 
  geom_point(data = random_distributions_1000, aes(X97.72, size), color = "palegreen3") + 
  geom_point(data = random_distributions_1000, aes(X84.13, size), color = "steelblue2") +
  geom_point(data = random_distributions_1000, aes(X99.86, size), color = "sienna2") +
  geom_point() + scale_color_manual(values = c("darkgrey", "blue")) +
  theme_light(base_size = 14)
plot(random_genes_plot_coloursig)

# save image of the R session
projectName <- "abortiveinfection_TPMamend"
save.image(file=paste0(projectName, ".RData"))
