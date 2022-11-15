library(phyloseq)
library(dplyr)
library(stringr)
library(magrittr)
library(data.table)
library(ggplot2)
library(ggsci)

###Data Import and Pre-processing###

OTU.matrix <- data.matrix(read.table(file="Sterilization_OTU.shared", header=TRUE, row.names=1, sep="\t"))
tax.matrix <- as.matrix(read.table(file="Sterilization_OTU.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="Sterilization.meta.csv", header=TRUE, row.names=1, sep=",")

#make joined OTU_Genus names
tax.data <- data.table(tax.matrix, keep.rownames = "Otu")
tax.data$Species2 <- paste(tax.data$Otu, tax.data$Genus, sep = "_")
tax.matrix <- as.matrix(data.frame(tax.data, row.names = "Otu"))

#Make nicer names for figures
sample <- sample %>% mutate (across('RawPulse', str_replace, 'chickpeaflour', 'Chickpea')) %>% 
  mutate (across('RawPulse', str_replace, 'sproutedgreenlentilflour', 'Sprouted_Green_Lentil')) %>%
  mutate (across('RawPulse', str_replace, 'greenlentilflour', 'Green_Lentil')) %>%
  mutate (across('RawPulse', str_replace, 'fieldpeaflour', 'Field_Pea'))

#Create phyloseq object
OTU.p <- otu_table(OTU.matrix, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
physeq <- phyloseq(OTU.p,TAX.p,sample.p)

#Remove irrelevant samples
physeq <- subset_samples(physeq, Media == "1XRUM")
physeq <- subset_samples(physeq, RawPulse != "nopulse")
physeq <- subset_samples(physeq, Dish != "15mlfalcontube")


##Endogenous microbe fermentations (no fecal samples), unsterilized pulses only

physeq_nf <- subset_samples(physeq, Fecal == "nofecal")
physeq_nf_un <- subset_samples(physeq_nf, Sterilization == "nosterilization")

#narrow down to taxa that are at least 0.1% abundant
physeq_nf_hi <- filter_taxa(physeq_nf_un, function(x) sum(x) > (sum(taxa_sums(physeq_nf_un))/1000), prune=TRUE)
#transform to relative abundances and account for the 2 fermentations per condition (divide by 2)
physeq_nf_rel <- transform_sample_counts(physeq_nf_hi, function(x) x/sum(x)/2)

tiff("Species_nf_bar_point1percent_unsterilized.tiff", units="in", width=8, height=5, res=300)
plot_bar(physeq_nf_rel, x= "RawPulse", fill="Species2") +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_d3("category20c", name = "Taxa") +
  xlab("Pulse Flour") +
  ylab("Relative Abundance")
dev.off()


##Simulated colonic fermentations (include fecal samples)
physeq_f <- subset_samples(physeq, Fecal == "M_001")

physeq_f_phyl <- tax_glom(physeq_f, taxrank = "Phylum")
physeq_f_fam <- tax_glom(physeq_f, taxrank = "Family")
physeq_f_genus <- tax_glom(physeq_f, taxrank = "Genus")

#narrow down to most abundant taxa
physeq_f_fam_hi <- filter_taxa(physeq_f_fam, function(x) sum(x) > (sum(taxa_sums(physeq_f_fam))/100), prune=TRUE)
physeq_f_gen_hi <- filter_taxa(physeq_f_genus, function(x) sum(x) > (sum(taxa_sums(physeq_f_genus))/100), prune=TRUE)

#Check mid-level taxa
physeq_f_fam_med <- filter_taxa(physeq_f_fam, function(x) sum(x) > (sum(taxa_sums(physeq_f_fam))/1000), prune=TRUE)
physeq_f_fam_med <- filter_taxa(physeq_f_fam_med, function(x) sum(x) < (sum(taxa_sums(physeq_f_fam))/100), prune=TRUE)

physeq_f_gen_med <- filter_taxa(physeq_f_genus, function(x) sum(x) > (sum(taxa_sums(physeq_f_genus))/1000), prune=TRUE)
physeq_f_gen_med <- filter_taxa(physeq_f_gen_med, function(x) sum(x) < (sum(taxa_sums(physeq_f_genus))/100), prune=TRUE)

#transform to relative abundances and account for the 2 fermentations per condition (divide by 2)
physeq_f_phyl_rel <- transform_sample_counts(physeq_f_phyl, function(x) x/sum(x)/2)

physeq_f_fam_rel <- transform_sample_counts(physeq_f_fam_hi, function(x) x/sum(x)/2)
physeq_f_fam_med_rel <- transform_sample_counts(physeq_f_fam_med, function(x) x/sum(x)/2)

physeq_f_gen_rel <- transform_sample_counts(physeq_f_gen_hi, function(x) x/sum(x)/2)
physeq_f_gen_med_rel <- transform_sample_counts(physeq_f_gen_med, function(x) x/sum(x)/2)

tiff("Phylum_f_bar.tiff", units="in", width=8, height=5, res=300)
plot_bar(physeq_f_phyl_rel, x= "Sterilization", fill="Phylum", facet_grid =~RawPulse) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_d3("category20c")
dev.off()

tiff("Family_f_bar_1percent.tiff", units="in", width=8, height=5, res=300)
plot_bar(physeq_f_fam_rel, x= "Sterilization", fill="Family", facet_grid =~RawPulse) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_d3("category20c")
dev.off()

tiff("Family_f_bar_point1to1.tiff", units="in", width=8, height=5, res=300)
plot_bar(physeq_f_fam_med_rel, x= "Sterilization", fill="Family", facet_grid =~RawPulse) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_d3("category20c")
dev.off()

tiff("Genus_f_bar_1percent.tiff", units="in", width=8, height=5, res=300)
plot_bar(physeq_f_gen_rel, x= "Sterilization", fill="Genus", facet_grid =~RawPulse) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_d3("category20c")
dev.off()

tiff("Genus_f_bar_poin1to1percent.tiff", units="in", width=8, height=5, res=300)
plot_bar(physeq_f_gen_med_rel, x= "Sterilization", fill="Genus", facet_grid =~RawPulse) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_d3("category20c")
dev.off()