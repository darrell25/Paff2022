library(data.table)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(multcomp)
library(DescTools)
library(phyloseq)
library(pairwiseAdonis)
library(DivNet)
library(vegan)
library(microbiome)
library(magrittr)
library(stringr)

###Data Import and Pre-processing###

OTU.matrix <- data.matrix(read.table(file="Sterilization_OTU.shared", header=TRUE, row.names=1, sep="\t"))
tax.matrix <- as.matrix(read.table(file="Sterilization_OTU.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="Sterilization.meta.csv", header=TRUE, row.names=1, sep=",")

#Create phyloseq object
OTU.p <- otu_table(OTU.matrix, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
tree.p <- read_tree(treefile="Sterilization.tree")
physeq <- phyloseq(OTU.p,TAX.p,sample.p, tree.p)

#Condense to genus level
physeq_gen <- tax_glom(physeq, taxrank="Genus")

#subset data to get appropriate samples for analysis
physeq_subset <- subset_samples(physeq_gen, Media=="1XRUM")
physeq_subset <- subset_samples(physeq_subset, Fecal=="M_001")
physeq_subset <- subset_samples(physeq_subset, RawPulse != "nopulse")
physeq_subset <- subset_samples(physeq_subset, Dish == "96deepwellplate")

###Alpha Diversity###
#Utilize the divnet package to compute alpha diversity metrics
a_div <- physeq_subset %>% divnet
sample.means <- data.frame(sample_data(physeq_subset))
sample.means$Shannon <- a_div$shannon %>% summary %$% estimate
sample.means$Simpson <- a_div$simpson %>% summary %$% estimate
sample.means$Simpson <- 1/sample.means$Simpson

#Modify names within sterilization column to improve readability of axis in graphic outputs:
sample.means <- sample.means %>% mutate (across('Sterilization', str_replace, "1.0bleach", "1.0% Bleach")) %>% 
  mutate (across('Sterilization', str_replace, "1.5bleach", "1.5% Bleach")) %>%
  mutate (across('Sterilization', str_replace, "70reagentalcohol", "70% Alcohol")) %>%
  mutate (across('Sterilization', str_replace, "33minuteautoclave", "Autoclave")) %>%
  mutate (across('Sterilization', str_replace, "autoclave", "Autoclave")) %>%
  mutate (across('Sterilization', str_replace, "2.0H2O2", "2.0% H2O2")) %>%
  mutate (across('Sterilization', str_replace, "nosterilization", "No Sterilization"))

#Convert sterilization values to factors and set No Sterilization as reference
Sterilizations <- c("Autoclave", "70% Alcohol", "1.5% Bleach", "1.0% Bleach", "2.0% H2O2", "No Sterilization", "none")

sample.means$Sterilization <- factor(sample.means$Sterilization, levels = Sterilizations, labels = Sterilizations)
sample.means$Sterilization <- relevel(as.factor(sample.means$Sterilization), ref = "No Sterilization")

#Linear model of Shannon diversity vs Sterilization with Dunnett test and BH correction
LinMod <- lm(Shannon ~ Sterilization, data = sample.means)
but_T <- glht(LinMod, linfct = mcp(Sterilization = "Dunnett"))
summary(but_T, test=adjusted("BH"))

#Linear model of Inverse Simpson diversity vs Sterilization with Dunnett test and BH correction
LinMod <- lm(Simpson ~ Sterilization, data = sample.means)
but_T <- glht(LinMod, linfct = mcp(Sterilization = "Dunnett"))
summary(but_T, test=adjusted("BH"))

#Generate graphics
tiff("Shannon_plot.tiff", units="in", width=12, height=5, res=300)
ggboxplot(sample.means, x = "Sterilization", y = "Shannon", add = "point")+
  xlab("Sterilization Treatment")+
  scale_x_discrete(limits = c("70% Alcohol","Autoclave","1.0% Bleach", "1.5% Bleach", "2.0% H2O2", "No Sterilization"))
dev.off()

tiff("Simpson_plot.tiff", units="in", width=12, height=5, res=300)
ggboxplot(sample.means, x = "Sterilization", y = "Simpson", add = "point")+
  ylab("Inverse Simpson")+
  xlab("Sterilization Treatment")+
  scale_x_discrete(limits = c("70% Alcohol","Autoclave","1.0% Bleach", "1.5% Bleach", "2.0% H2O2", "No Sterilization"))
dev.off()

###Beta-Diversity###

##Bray-Curtis##

#Calculate Bray-dissimilarity
otu.mat <- t(as(otu_table(physeq_subset), "matrix"))
b_veg_bray <- vegdist(otu.mat, method = "bray")
b_veg_bray_pcoa <- ordinate(physeq_subset, method="PCoA", distance=b_veg_bray)

#Perform PERMANOVA
set.seed(3488)
adonis2(formula = b_veg_bray ~ Sterilization, data = sample.means)

#Perform Pair-wise PERMANOVA to find differences from unsterilized
set.seed(3488)
pairwise.adonis2(b_veg_bray ~ Sterilization, data = sample.means, nperm = 999)

tiff("Bray_plot.tiff", units="in", width=5, height=5, res=300)
phyloseq::plot_ordination(physeq_gen, b_veg_bray_pcoa, type ="samples", color = "Sterilization") +
  ggtitle("Bray-Curtis Dissimilarity") +
  theme (plot.title = element_text(hjust = 0.5))+
  scale_color_discrete(name = "Treatment", labels = c("1.0% Bleach", "1.5% Bleach", "2.0% Hydrogen Peroxide", 
                                                      "Autoclave", "70% Reagent Alcohol", "No Sterilization"))
dev.off()

##Aitchison##

#CLR transformation of counts
physeq_clr <- microbiome::transform(physeq_subset, "clr")
otu.clr <- t(as(otu_table(physeq_clr), "matrix"))

#calculate aitchison distance (Euclidan distance of clr transformed data)
clr_ait <- dist(otu.clr, method='euc')
clr_ait_pcoa <- ordinate(physeq_clr, method="PCoA", distance=clr_ait)

#Perform PERMANOVA
set.seed(3488)
adonis2(formula = clr_ait ~ Sterilization, data = sample.means)

#Perform Pair-wise PERMANOVA to find differences from unsterilized
set.seed(3488)
pairwise.adonis2(clr_ait ~ Sterilization, data = sample.means, nperm = 999)

#Generate graphic
tiff("Aitchison_Plot.tiff", units="in", width=5, height=5, res=300)
phyloseq::plot_ordination(physeq_clr, clr_ait_pcoa, type ="samples", color = "Sterilization") +
  ggtitle("Aitchison Distance") +
  theme (plot.title = element_text(hjust = 0.5))+
  scale_color_discrete(name = "Treatment", labels = c("1.0% Bleach", "1.5% Bleach", "2.0% Hydrogen Peroxide", 
                                                      "Autoclave", "70% Reagent Alcohol", "No Sterilization"))
dev.off()

##Weighted UniFrac##

#Rarify the dataset (recommended for UniFrac)
set.seed(3488)
physeq_subset.rar <- rarefy_even_depth(physeq_subset, replace = FALSE)

#Calculate UniFrac Beta Diversity:
dist_gen_rar_uniw <- phyloseq::UniFrac(physeq_subset.rar, weighted = TRUE) 
gen_rar_uniw_pcoa <- ordinate(physeq_subset.rar, method="PCoA", distance=dist_gen_rar_uniw)

#Perform PERMANOVA
set.seed(3488)
adonis2(formula = dist_gen_rar_uniw ~ Sterilization, data = sample.means)

#Perform Pair-wise PERMANOVA to find differences from unsterilized
set.seed(3488)
pairwise.adonis2(dist_gen_rar_uniw ~ Sterilization, data = sample.means, nperm = 999)


tiff("UniFrac_Plot.tiff", units="in", width=5, height=5, res=300)
phyloseq::plot_ordination(physeq_subset.rar, gen_rar_uniw_pcoa, type ="samples", color = "Sterilization") +
  ggtitle("Weighted UniFrac") +
  theme (plot.title = element_text(hjust = 0.5))+
  scale_color_discrete(name = "Treatment", labels = c("1.0% Bleach", "1.5% Bleach", "2.0% Hydrogen Peroxide", 
                                                      "Autoclave", "70% Reagent Alcohol", "No Sterilization"))
dev.off()
