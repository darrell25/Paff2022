library(DESeq2)
library(phyloseq)
library(data.table)

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

###DESeq2 Analysis###

#convert to factors
sample2 <- data.frame(sample_data(physeq_83_5))
sample2$Sterilization <- factor(sample2$Sterilization, levels = 
                                  c("nosterilization","33minuteautoclave", 
                                    "70reagentalcohol", "2.0H2O2", "1.0bleach", "1.5bleach" ))

sterile <- factor(c("33minuteautoclave", "70reagentalcohol", "2.0H2O2", "1.0bleach", "1.5bleach"))

#Run analysis looking at differences in Sterilization
OTU83_5.matrix <- as.data.frame(otu_table(physeq_83_5))
DE_Diet_data83_5 <- DESeqDataSetFromMatrix(countData = OTU83_5.matrix, colData = sample2, 
                                           design = ~Sterilization)
DE_Diet83_5 <- DESeq(DE_Diet_data83_5)

#find the results comparing sterilization treatment to unsterilized
for (fac in sterile){
  res83_5 <- results(DE_Diet83_5, contrast = c("Sterilization", "nosterilization", fac), cooksCutoff = FALSE)
  alpha_cut <- 0.05
  sigtab83_5 <- res83_5[which(res83_5$padj < alpha_cut), ]
  sigtab83_5 = cbind(as(sigtab83_5, "data.frame"), as(tax_table(physeq_83_5)[rownames(sigtab83_5), ], "matrix"))
  
  fn <- paste("DESEQ_", fac, ".csv", sep = "")
  write.csv(sigtab83_5,file=fn)
}
