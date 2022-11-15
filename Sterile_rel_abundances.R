library(phyloseq)
library(data.table)
library(dplyr)

##Make Phyloseq Object##

#Import data
OTU.matrix <- data.matrix(read.table(file="Sterilization_OTU.shared", header=TRUE, row.names=1, sep="\t"))
tax.dt <- fread(file = "Sterilization_OTU.taxonomy", header = TRUE, sep="\t")
sample <- read.csv(file="Sterilization.meta.csv", header=TRUE, row.names=1, sep=",")

#Create joined genus_species names
tax.dt$Species <- make.names (paste(tax.dt$OTU, tax.dt$Genus, sep = "_"), unique=TRUE)

tax.matrix<-as.matrix(tax.dt, rownames = "OTU")

#Create phyloseq object
OTU.p <- otu_table(OTU.matrix, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
physeq <- phyloseq(OTU.p,TAX.p,sample.p)

##Filter sequences##

#Subset Data
physeq <- subset_samples(physeq, Media == "1XRUM")
physeq <- subset_samples(physeq, RawPulse != "nopulse")
physeq <- subset_samples(physeq, Fecal == "M_001")
physeq <- subset_samples(physeq, Dish != "15mlfalcontube")

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

#Get the average relative abundance for each of the diet conditions for each member 
physeq_83_5 <- transform_sample_counts(physeq_83_5, function(x) x/sum(x))

physeq_us <- data.frame (otu_table(subset_samples(physeq_83_5, Sterilization == "nosterilization")))
physeq_10B <- data.frame (otu_table(subset_samples(physeq_83_5, Sterilization == "1.0bleach")))
physeq_15B <- data.frame (otu_table(subset_samples(physeq_83_5, Sterilization == "1.5bleach")))
physeq_H2O2 <- data.frame (otu_table(subset_samples(physeq_83_5, Sterilization == "2.0H2O2")))
physeq_alc <- data.frame (otu_table(subset_samples(physeq_83_5, Sterilization == "70reagentalcohol")))
physeq_auto <- data.frame (otu_table(subset_samples(physeq_83_5, Sterilization == "33minuteautoclave")))

Sterile.pct <- as.data.table(data.frame(OTU=rownames(tax_table(physeq_83_5)), Species=tax_table(physeq_83_5)[,"Species"], 
                          Unsterilized = round(rowMeans(physeq_us)*100,3), 
                          Alcohol = round(rowMeans(physeq_alc)*100,3), Autoclave = round(rowMeans(physeq_auto)*100,3), 
                          Bleach_1.0 = round(rowMeans(physeq_10B)*100,3), Bleach_1.5 = round(rowMeans(physeq_15B)*100,3),
                          H2O2 = round(rowMeans(physeq_H2O2)*100,3)))

fwrite(Sterile.pct, file="OTU_percent_abundances.txt", sep = "\t")

