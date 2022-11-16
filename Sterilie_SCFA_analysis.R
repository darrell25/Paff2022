library(data.table)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(multcomp)
library(lme4)
library(lmerTest)
library(DescTools)
library(stringr)

Sterilizations <- c("Autoclave", "70% Alcohol", "1.5% Bleach", "1.0% Bleach", "2.0% H2O2", "No Sterilization", "none")

sample.means <- fread(file="sterilization.meta.csv")

#Modify names within sterilization column to improve readability of axis in graphic outputs:
sample.means <- sample.means %>% mutate (across('Sterilization', str_replace, "1.0bleach", "1.0% Bleach")) %>% 
  mutate (across('Sterilization', str_replace, "1.5bleach", "1.5% Bleach")) %>%
  mutate (across('Sterilization', str_replace, "70reagentalcohol", "70% Alcohol")) %>%
  mutate (across('Sterilization', str_replace, "33minuteautoclave", "Autoclave")) %>%
  mutate (across('Sterilization', str_replace, "autoclave", "Autoclave")) %>%
  mutate (across('Sterilization', str_replace, "2.0H2O2", "2.0% H2O2")) %>%
  mutate (across('Sterilization', str_replace, "nosterilization", "No Sterilization"))


sample.means$Sterilization <- factor(sample.means$Sterilization, levels = Sterilizations, labels = Sterilizations)
sample.means$Sterilization <- relevel(as.factor(sample.means$Sterilization), ref = "No Sterilization")

#Sub-setting what data gets included in analysis: 
sample.subset <- sample.means[Fecal == "M_001", ]
sample.subset <- sample.subset[Media == "1XRUM", ]
sample.subset <- sample.subset[RawPulse != "nopulse", ]
sample.subset <- sample.subset[Dish == "96deepwellplate"]

##Butyrate##
#Dunnett test with BH pairwise correction (below), looking at pairwise comparisons against not treating.
LinMod <- lm(Butyrate_mM ~ Sterilization, data = sample.subset)
but_T <- glht(LinMod, linfct = mcp(Sterilization = "Dunnett"))
summary(but_T, test=adjusted("BH"))

#Generate graphic
tiff("Butyrate.tiff", units="in", width=12, height=5, res=300)
ggboxplot(sample.subset, x = "Sterilization", y = "Butyrate_mM", add = "point")+
  xlab("Sterilization Treatment")+
  scale_x_discrete(limits = c("70% Alcohol","Autoclave","1.0% Bleach", "1.5% Bleach", "2.0% H2O2", "No Sterilization"))
dev.off()

##Acetate##
#Dunnett test with BH pairwise correction (below), looking at pairwise comparisons against not treating.
LinMod <- lm(Acetate_mM ~ Sterilization, data = sample.subset)
but_T <- glht(LinMod, linfct = mcp(Sterilization = "Dunnett"))
summary(but_T, test=adjusted("BH"))

#Generate graphic
tiff("Acetate.tiff", units="in", width=12, height=5, res=300)
  ggboxplot(sample.subset, x = "Sterilization", y = "Acetate_mM", add = "point")+
  xlab("Sterilization Treatment")+
  scale_x_discrete(limits = c("70% Alcohol","Autoclave","1.0% Bleach", "1.5% Bleach", "2.0% H2O2", "No Sterilization"))
dev.off()
