setwd("~/Documents/R/hepatitis_analysis/data")
###### This is modification of the original correlation script by Cristina Venturini

###### Liver transcriptome - hepatitis
###### 30/09/2022

library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)

########## 1. Average expression cluster target gene  ########## 
## A. Input files - replace as necessary

transcriptome<- read.csv("example_tpm_PC0.001_log2_genesymbol_dedup.csv", header=T, row.names = 1)
ipa <- read.delim("SigCorr_dedupURA_Z2_d2tst_dge-upAll_mwtest_fdr05.txt",header=T) # IPA output file

## B. Create a mean expression for each gene across individuals 
transcriptome <- transcriptome[order(rownames(transcriptome)),]
mean_df <- as.data.frame(rowMeans(transcriptome, na.rm=TRUE)) #mean by gene
colnames(mean_df) <- c("AverageSamples")

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
ipa_sep_reduced <- subset(ipa_sep, target %in% mean_df$target) #only the genes found in the transcriptome co-correlation matrix
#NB. sometimes on my computer this subsetting is off. If dim(ipa_sep_reduced) shows 0, try unload dplyr and load again. 
dim(ipa_sep_reduced)

ipa_sep_reduced %>%
  group_by(ID) %>%
  mutate("count" = n()) -> ipa_sep_reduced # add a column with the cluster size

## D. Calculate the average of expression for each cluster 

mean_df$target <- rownames(mean_df)
ipa_sep_reduced_avg <- left_join(ipa_sep_reduced,mean_df,by="target") #merge cluster df with avg by gene

ipa_avg_bycluster <- ipa_sep_reduced_avg %>% group_by(regulator,count,ID,type,FDR) %>% 
  summarise(AvgCluster = mean(AverageSamples))
#WARNING:`summarise()` has grouped output by 'regulator', 'count', 'ID', 'type'. You can override using the `.groups` argument. IT IS FINE. 

#ipa_avg_bycluster will include all info about the clusters 

dim(ipa_avg_bycluster)
#195 clusters 
## E. Output file
setwd("~/Documents/R/hepatitis_analysis/output")
write.csv(ipa_avg_bycluster, "1_avg-exp_hepatitis.csv", row.names = F)
###=============================================================================
########## 2. Select random clusters of genes  ########## 
# generate random distributions

possiblesizes <- unique(ipa_avg_bycluster$count) #clusters' size present in IPA

#function to pick genes randomly to create clusters the same size of IPA's clusters. 
#data_df is the dataframe containing all average expression values for each gene, N is the number of iterations (100), SIZECLUSTER a vector with the possible clusters' size (can be 4:807 or "possiblesizes" vector)
#map_dfr is a handy R function to avoid for loops (which can be slow in R)
samplingclusters <- function(data_df,N,SIZECLUSTER){
  map_dfr(seq_len(N), ~ data_df %>% sample_n(SIZECLUSTER) %>% 
            mutate(iteration = .x))  %>% 
    mutate(clustersize = SIZECLUSTER) %>% 
    group_by(iteration,clustersize) %>% 
    summarise(Avg = mean(AverageSamples))
}
#running the function - with 100 iteration this takes < 1 min
res <- possiblesizes %>%
  map_dfr(samplingclusters, data=mean_df,N=100)

#calculate mean for each clusters' size, sd and other 
random_distributions_100 <- res %>% group_by(clustersize) %>% 
  summarise(AvgExp=mean(Avg),SDExp=sd(Avg),X84.13=AvgExp+SDExp,X97.72=AvgExp+2*SDExp,X99.86=AvgExp+3*SDExp) 

#SOME PLOTS

ggplot(random_distributions_100, aes(x=AvgExp)) + 
  geom_histogram(fill = "red", alpha = 0.2) +
  geom_histogram(data=ipa_avg_bycluster,aes(x=AvgCluster),fill = "blue", alpha = 0.2) +
  #xlim(-9.965784,14.957397) +
  theme_classic()

ggplot(random_distributions_100,aes(x=AvgExp,y=clustersize)) +
  geom_point() +
  theme_classic()

#ipa_avg_bycluster
ggplot(ipa_avg_bycluster,aes(x=AvgCluster,y=count)) +
  geom_point() +
  theme_classic()


## D. Output file
#write.csv(random_distributions_100,
colnames(random_distributions_100) <- c("size", "mean", "sd", "X84.13", "X97.72",	"X99.86")
write.csv(random_distributions_100,
          paste0("2_random_distributions_hepatitis",Sys.Date(),".csv"),
          row.names = F)

#===================================================================
########## 3. Calc stats

ipa_avg_expression <- read.csv("1_avg-exp_hepatitis.csv") # Output file from Script 1

random_distributions_100 <- read.csv("2_random_distributions_hepatitis2022-10-03.csv") # output file from Script 2 

## B. Determining which clusters are FDR significant
ipa_avg_expression <- ipa_avg_expression[ipa_avg_expression$count >= 4,] # take out clusters smaller than 4
ipa_avg_expression$rc_mean <- random_distributions_100$mean[match(ipa_avg_expression$count, random_distributions_100$size)] #mean
ipa_avg_expression$rc_sd <- random_distributions_100$sd[match(ipa_avg_expression$count, random_distributions_100$size)] #sd

ipa_avg_expression$zscore <- (ipa_avg_expression$AvgCluster - ipa_avg_expression$rc_mean)/ipa_avg_expression$rc_sd #zscores

ipa_avg_expression$pvalue <- pnorm(ipa_avg_expression$zscore, lower.tail = FALSE) #p-values
ipa_less_0.05 <- ipa_avg_expression[which(ipa_avg_expression$pvalue <= 0.05),] # p-values <= 0.05
fdr_threshold <- nrow(ipa_less_0.05)
ipa_less_0.05$adj_pvalue <- ipa_less_0.05$pvalue*fdr_threshold # calculating the adjusted p-value
ipa_significant <- ipa_less_0.05[which(ipa_less_0.05$adj_pvalue <= 0.05),] # selecting clusters with adj. p-value <= 0.05

## C. Output files
#write.csv(ipa_significant, "../data/UR_module_identification/3_fdrsig_40size_inf_pre0.csv", row.names = F) # clusters with adj. p-value <= 0.05
write.csv(ipa_significant, "3_significantexpression_hepatitisB.csv", row.names = F) # clusters with adj. p-value <= 0.05
write.csv(ipa_less_0.05, "3_expression_ipaless0_05_hepatitisB.csv", row.names = F)
write.csv(ipa_avg_expression, "3_ipa_avgexpression_hepatitisB.csv", row.names = F)
## D. Plotting the FDR-significant clusters
library(ggplot2)

# plotting all the clusters against the frequency distribution and making FDR-sig clusters blue
ipa_avg_expression$Significance <- (ipa_avg_expression$regulator%in%ipa_significant$regulator) # need to have a T/F significance column
random_genes_plot_coloursig <- ggplot(ipa_avg_expression, aes(AvgCluster, count, colour = Significance)) + 
  xlab("Average expression") + ylab("Network size (number of target genes)") + 
  geom_point(data = random_distributions_100, aes(X97.72, size), color = "palegreen3") + 
  geom_point(data = random_distributions_100, aes(X84.13, size), color = "steelblue2") +
  geom_point(data = random_distributions_100, aes(X99.86, size), color = "sienna2") +
  geom_point() + scale_color_manual(values = c("darkgrey", "blue")) +
  theme_light(base_size = 14)
plot(random_genes_plot_coloursig)

# save image of the R session
projectName <- "AVG_EXPRESSION_hepatitis"
save.image(file=paste0(projectName, ".RData"))


#do again for all clust zscores

ipa_avg_expression <- read.csv("1_avg-exp_hepatitis.csv") # Output file from Script 1

random_distributions_100 <- read.csv("2_random_distributions_hepatitis2022-10-03.csv") # output file from Script 2 

## B. Determining which clusters are FDR significant
ipa_avg_expression$rc_mean <- random_distributions_100$mean[match(ipa_avg_expression$count, random_distributions_100$size)] #mean
ipa_avg_expression$rc_sd <- random_distributions_100$sd[match(ipa_avg_expression$count, random_distributions_100$size)] #sd

ipa_avg_expression$zscore <- (ipa_avg_expression$AvgCluster - ipa_avg_expression$rc_mean)/ipa_avg_expression$rc_sd #zscores

ipa_avg_expression$pvalue <- pnorm(ipa_avg_expression$zscore, lower.tail = FALSE) #p-values
fdr_threshold <- nrow(ipa_avg_expression)
ipa_avg_expression$adj_pvalue <- ipa_avg_expression$pvalue*fdr_threshold # calculating the adjusted p-value


## C. Output files
#write.csv(ipa_significant, "../data/UR_module_identification/3_fdrsig_40size_inf_pre0.csv", row.names = F) # clusters with adj. p-value <= 0.05

write.csv(ipa_avg_expression, "3_ipa_avgexpression_nomliv_allmod.csv", row.names = F)
## D. Plotting the FDR-significant clusters
