library(exactRankTests)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library(phyloseq)
library(data.table)
library(tidyverse)
library(dplyr)
library(magrittr)

# ANCOM 2.1 from: https://github.com/FrederickHuangLin/ANCOM/blob/master/scripts/ancom_v2.1.R
# Place in working directory
source("ancom_v2.1.R")

###Data Import and Processing###
OTU.matrix <- data.matrix(read.table(file="Sterilization_OTU.shared", header=TRUE, row.names=1, sep="\t"))
tax.matrix <- as.matrix(read.table(file="Sterilization_OTU.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="Sterilization.meta.csv", header=TRUE, row.names=1, sep=",")

#make joined OTU_Genus names
tax.data <- data.table(tax.matrix, keep.rownames = "Otu")
tax.data$Species2 <- paste(tax.data$Otu, tax.data$Genus, sep = "_")
tax.matrix <- as.matrix(data.frame(tax.data, row.names = "Otu"))

##Create phyloseq object##
OTU.p <- otu_table(OTU.matrix, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
physeq <- phyloseq(OTU.p,TAX.p,sample.p)

##Subset Data##
#Remove samples that are not relevant for this analysis
physeq <- subset_samples(physeq, Media == "1XRUM")
physeq <- subset_samples(physeq, RawPulse != "nopulse")
physeq <- subset_samples(physeq, Fecal == "M_001")
physeq <- subset_samples(physeq, Dish != "15mlfalcontube")

##Filter sequences##

#Create filter to have minimum 0.001% of total reads and 5% prevalence

physeq_83 <- filter_taxa(physeq, function(x) sum(x) >=(sum(taxa_sums(physeq))*0.00001), prune=TRUE)

# Compute prevalence of each feature, store as data.frame
prevdf83 <- apply(X = otu_table(physeq_83),
                  MARGIN = ifelse(taxa_are_rows(physeq_83), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf83 <- data.frame(Prevalence = prevdf83,
                       TotalAbundance = taxa_sums(physeq_83),
                       tax_table(physeq_83))

prevalenceThreshold <- 0.05 * nsamples(physeq_83)
keepTaxa2 <- rownames(prevdf83)[(prevdf83$Prevalence >= prevalenceThreshold)]
physeq_83_5 <- prune_taxa(keepTaxa2, physeq_83)


#setup variables for ANCOM
sample_var <- c("SeqNum")
lib_cut = 1000
neg_lb=FALSE   
main_var <- c("Sterilization")
sterile <- factor(c("33minuteautoclave", "70reagentalcohol", "2.0H2O2", "1.0bleach", "1.5bleach"))

for (fac in sterile){
  #set pair for comparison, change here and output file name
  physeq_test <- subset_samples(physeq_83_5, Sterilization %in% c("nosterilization", fac))
  
  #Data Pre-Processing
  feature_table <- as.data.frame(otu_table(physeq_test))
  meta_data <- as.data.table(data.frame(sample_data(physeq_test)), keep.rownames = "SeqNum")
  processed_data <- feature_table_pre_process(feature_table, meta_data, sample_var, group_var=NULL, 
                                              out_cut = 0.05, zero_cut = 0.95, lib_cut, neg_lb)
  
  #Run ANCOM
  ANCOM_results <- ANCOM(processed_data$feature_table, processed_data$meta_data, 
                         struc_zero = NULL, main_var, p_adj_method = "BH", 
                         alpha = 0.05, adj_formula = NULL, rand_formula = NULL)
  
  #Process results to get only the significantly different OTUs and reformat for output
  results_83_5 <- data.table(ANCOM_results$out)
  results_0.6 <- results_83_5[detected_0.6==TRUE,]
  result_dat <- ANCOM_results$fig$data
  
  theresults <- filter(result_dat, taxa_id %in% results_0.6$taxa_id)
  theresults <- cbind(theresults, results_0.6[,detected_0.9:detected_0.6])
  theresults <- dplyr::rename(theresults, meanCLR=x)
  theresults <- dplyr::rename(theresults, W_score=y)
  
  spec_nam <- as.data.table(as.data.frame(tax_table(physeq_test)), keep.rownames = "OTU")
  theresults <- cbind(taxa_id=theresults$taxa_id, Species=spec_nam[spec_nam$OTU %in% theresults$taxa_id, Species], 
                      theresults[,c(2:8)])
  
  fn <- paste("ANCOM_", fac, ".sig", sep = "")
  fwrite(theresults, file = fn, sep = "\t") 
}

