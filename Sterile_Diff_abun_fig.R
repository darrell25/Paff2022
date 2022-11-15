library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)

#Import data
da_tab <- fread("Diff_abund.csv")
rel_sp <- fread("OTU_percent_abundances.txt")

#Get all of the significant relative abundances
sp_names <- da_tab$OTU
rel_sp_sig <- rel_sp[OTU %in% sp_names,]

#Reformat to facilitate reusing some code
setnames(rel_sp_sig, "Species", "Name")
rel_sp_sig <- rel_sp_sig[,-1]
rel_all <- rel_sp_sig

#Calculate normalized values for differences between sterilized and unsterilized
rel_nums <- rel_all[,-c(1:2)]
unsterilized <- rel_all$Unsterilized
diff_nums <- rel_nums - unsterilized
norm_nums <- diff_nums/apply(abs(diff_nums), 1, max)

#Prepare objects for plotting
plot_vals <- cbind(Name=rel_all$Name, norm_nums)
plot_sigs <- cbind(Name=rel_all$Name, da_tab[,c(3:7)])
colnames(plot_sigs) <- colnames(plot_vals) #don't think its needed

#set the order of names by OTU number (have to do this before melting data to avoid repeat names)
ordernames <- sort(plot_vals$Name, decreasing = TRUE) 

#convert data to proper format for plotting (long)
plot_vals <- data.table::melt(plot_vals, variable.name = "Sterilization", id.vars = "Name", value.name = "Rel_Change")
plot_sigs <- data.table::melt(plot_sigs, variable.name = "Sterilization", id.vars = "Name", value.name = "Sig")

#Convert sterilizations to factors to give me control over order
treat_ord <- c("Alcohol", "Autoclave", "Bleach_1.0", "Bleach_1.5", "H2O2")
plot_vals$Sterilization <- factor(plot_vals$Sterilization, levels = treat_ord)
plot_sigs$Sterilization <- factor(plot_sigs$Sterilization, levels = treat_ord)

#Convert names to factors to give me control over order
plot_vals$Name <- factor(plot_vals$Name, levels = ordernames)
plot_sigs$Name <- factor(plot_sigs$Name, levels = ordernames)

#get the significant changes to add in the labels
plot_sigs <- plot_sigs[is.na(plot_sigs$Sig)==FALSE,]

sig_dat <- data.frame(Sterilization = plot_sigs$Sterilization, Name = plot_sigs$Name)

tiff("Diff_abundance2.tiff", units="in", width=7, height=6, res=300)
ggplot(plot_vals, aes(x = Sterilization, y = Name)) +
  geom_tile(aes(fill = Rel_Change)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", name = "Change From Control") +
  labs(y = "Taxa", x = "Sterilization Treatment") +
  scale_x_discrete(labels = c("Alcohol" = "70%_Alcohol", "Autoclave" = "Autoclave", "Bleach_1.0" = "1.0%_Bleach", 
                              "Bleach_1.5" = "1.5%_Bleach", "H2O2" = "2.0%_H2O2")) +
  geom_text(data = sig_dat, label = plot_sigs$Sig)
dev.off()
