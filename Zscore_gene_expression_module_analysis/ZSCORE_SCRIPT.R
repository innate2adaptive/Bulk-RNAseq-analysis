#####calculating z-scores 

#AC 2022

#set working directory first

#load necessary libraries 
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(data.table)


#load necessary files

#module file
ipa <- read.csv("celltypemodules.csv",header=T)
ipa <- tibble::rowid_to_column(ipa, "ID") #add an ID column corresponding to each cluster
colnames(ipa) <- c("ID","regulator","target") #rename columns

#load transcriptome file
transcriptome<- read.delim("GSE96851_AdultHepB.txt", header=T)
#remove duplicated genes
transcriptome <- transcriptome[!duplicated(transcriptome$Gene.Symbol),]
#rename first columns
colnames(transcriptome)[1] <- "target"


ipa$target <- as.character(ipa$target) #target column is integer -> convert to character for next step
ipa_sep <- separate_rows(ipa, target, sep = ",|/", convert = TRUE) #separate clusters into individual rows
ipa_sep_reduced <- subset(ipa_sep, target %in% transcriptome$target)

ipa_sep_reduced %>%
  group_by(ID) %>%
  mutate("count" = n()) -> ipa_sep_reduced # add a column with the cluster size


#join transcriptome and module file by target file 
ipa_sep_reduced_avg <- left_join(ipa_sep_reduced, transcriptome,by="target")

#check dimensions if joined correctly
dim(ipa_sep_reduced_avg)
head(ipa_sep_reduced_avg)


#calculate the average expression per gene (row) per column (sample)
means <- setDT(ipa_sep_reduced_avg)[, lapply(.SD, mean),.SDcols=colnames(ipa_sep_reduced_avg[,5:21]), by=regulator]

#again check the table 
dim(means)
means

#save cluster size into a separate vector and remove any repeated values
counts <- ipa_sep_reduced_avg[,2:4]
counts <- counts[,-2]
counts <- counts[!duplicated(counts$regulator),]
head(counts)

#join with the average expression table
means<- left_join(means, counts,by="regulator")
means


#set working directory and save into a file
write.csv(means,"meanexpr_persample_hep_unfil.csv")



#load your file with mean expression per sample calculated
ipa_avg_expression <- read.csv("meanexpr_persample_hep_unfil.csv") # Output file from Script 1

#load random distribution file for average expression with same cluster sizes
random_distributions_100 <- read.csv("exp_random_distributions_hepatitis_allmod_unfil2022-10-05.csv") # output file from Script 2 


## B. Determining which clusters are FDR significant

#take out average expression for 1st sample and save it into a vector
ipa_s1 <- ipa_avg_expression[,2:3]
#take cluster size and save into a vector
ipa_s1$count <- ipa_avg_expression$count

#match the cluster size and take mean expression from random distribution
ipa_s1$rc_mean <- random_distributions_100$mean[match(ipa_s1$count, random_distributions_100$size)] #mean
#extract standard deviation from random distribution file by matching for cluster size with sample
ipa_s1$rc_sd <- random_distributions_100$sd[match(ipa_s1$count, random_distributions_100$size)] #sd

#calculate z-score by taking (mean expression of sample 1 - mean expression of random same cluster size)/stdev of the random distribution
ipa_s1$zscore <- (ipa_s1$X22LM.117V0269.2_RNA_S1 - ipa_s1$rc_mean)/ipa_s1$rc_sd #zscores

#repeat the code for all samples
#due to time pressure, I've done this calculation manually across the entire data set, however I am currently working towards writing this as a function that could calculate this faster by looping through the dataset

#create new matrix to save results
zscores <- matrix(nrow=202,ncol=5)

#get names of regulators from the avg expression data set
zscores[,1] <- ipa_avg_expression$regulator
#add all zscores for each sample
zscores[,2] <- ipa_s1$zscore
zscores[,3] <- ipa_s2$zscore
zscores[,4] <- ipa_s3$zscore
zscores[,5] <- ipa_s4$zscore

##until you reach the last one - will be written as a function for automation as well
zscores[,18] <- ipa_s17$zscore

#first column is regulator 
colnames(zscores)[1] <- "regulator"
#copy sample names from a transcriptome file accordingly - to name columns
colnames(zscores)[2:18] <- colnames(transcriptome)[1:17]


#save onto a file
write.csv(zscores,"zscores_persample_hep_unfil.csv")
