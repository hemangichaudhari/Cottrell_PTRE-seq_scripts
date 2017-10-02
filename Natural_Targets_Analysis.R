library(dplyr)
library(ggplot2)
library(extrafont)
font_import()

loadfonts(device="win")

cbPalette <- c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_science <- function (base_size = 12, base_family = "Arial Black") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=2), 
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour= "black", size=1),  axis.line.y = element_line(colour= "black", size=1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}
#load RNA results
RNA_results= read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/RNA_plasmid.txt")
RNA_results_natural <- RNA_results
#load grouped_barcode identities
Grouped_barcode_identities = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Grouped_Barcode_identities.txt")
#filter RNA results with only those contianing Natural Let-7 binding sites
RNA_results_natural <- dplyr::filter(RNA_results_natural, RE_Identity %in% Grouped_barcode_identities$Natural)

#rename '7777' as 'Synthetic'
RNA_results_natural$RE_Identity <- gsub("7777", "Synthetic", RNA_results_natural$RE_Identity)

#plot RNA expression
ggplot(RNA_results_natural, aes(x = reorder(RE_Identity, -FOLD), FOLD)) + 
  xlab("") + ylab("Relative Expression (log2)") +  
  geom_boxplot(outlier.shape = 1, notch =TRUE) + theme_science() + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) 
ggsave("Natural_RNA.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#plot RNA expression with invisible outliers
natural_rna <- ggplot(RNA_results_natural, aes(x = reorder(RE_Identity, -FOLD), FOLD)) +  
  xlab("") + ylab("Relative Expression (log2)") + 
  geom_boxplot(outlier.colour = NA, notch =TRUE) + theme_science() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")
#rescale above plot
sts_natural_1 <- boxplot.stats(RNA_results_natural$FOLD)$stats
p1 = natural_rna + coord_cartesian(ylim = c(min(sts_natural_1), max(sts_natural_1)))
p1
ggsave("Natural_RNA_NO.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)


#new data.frame with MFE for base-pairing (MFE) and MFE for 3'UTR (MFE_3)
MFE <- data.frame(RE_Identity=c("Synthetic", "SMARCAD1", "C14orf28", "Dna2", "FIGNL2", "HMGA2"), 
                  MFE=c(-28.8, -20.0, -18.2, -16.9, -25.0, -23.5),
                  MFE_3=c(-87.7, -82.6, -72.1, -70.6, -82.3, -63.1))


#Merge RNA results for natural targets with MFE
Natural_MFE <-merge(RNA_results_natural, MFE, by="RE_Identity")
#drop control
Natural_MFE <- subset(Natural_MFE, !RE_Identity =="Control")
#group and summarise to get median Fold change of RNA
Natural_MFE <- dplyr::group_by(Natural_MFE, RE_Identity)
Natural_MFE <- dplyr:: summarise(Natural_MFE, Fold=median(FOLD), MAD=mad(FOLD), MFE=mean(MFE), MFE_3=mean(MFE_3))

#correlation between MFE and Fold
ggplot(Natural_MFE, aes(x=MFE, y=Fold)) +
  geom_point() 
cor.test(x= Natural_MFE$Fold, y= Natural_MFE$MFE, method="spearman")

# Formula for the linear model with MFE
mfe <- as.formula(paste("Fold ~ MFE"))
# Fit all data
fit_mfe <- lm(mfe, data = Natural_MFE)
summary(fit_mfe) #displays summary to the fit. Information stored in this variable
predpr<- predict(fit_mfe,type=c("response")) # Predicted values from the fit 
Natural_MFE$problm=predpr # Adding predicted values to the dataset
t= cor(Natural_MFE$problm, Natural_MFE$Fold)
ggplot(Natural_MFE, aes(Fold,problm)) + geom_point() + 
  ggtitle(label=(paste("R =", format(t, digits = 2)))) + 
  xlab("Measured Relative RNA Expression (log2)") + 
  ylab("Predicted Relative RNA Expression (log2)") + 
  geom_smooth(method = lm) +
  theme_science()

# Formula for the linear model with MFE 3'UTR
mfe_3 <- as.formula(paste("Fold ~ MFE_3"))
# Fit all data
fit_mfe_3 <- lm(mfe_3, data = Natural_MFE)
summary(fit_mfe_3) #displays summary to the fit. Information stored in this variable
predpr<- predict(fit_mfe_3,type=c("response")) # Predicted values from the fit 
Natural_MFE$problm_2=predpr # Adding predicted values to the dataset
t= cor(Natural_MFE$problm_2, Natural_MFE$Fold)
ggplot(Natural_MFE, aes(Fold,problm_2)) + geom_point() + 
  ggtitle(label=(paste("R =", format(t, digits = 2)))) + 
  xlab("Measured Relative RNA Expression (log2)") + 
  ylab("Predicted Relative RNA Expression (log2)") + 
  geom_smooth(method = lm) +
  theme_science()


# Formula for the linear model with both MFE and MFE 3'UTR
mfe_mfe_3 <- as.formula(paste("Fold ~ MFE*MFE_3"))
# Fit all data
fit_mfe_mfe_3 <- lm(mfe_mfe_3 , data = Natural_MFE)
summary(fit_mfe_mfe_3) #displays summary to the fit. Information stored in this variable
predpr<- predict(fit_mfe_mfe_3,type=c("response")) # Predicted values from the fit 
Natural_MFE$problm_3=predpr # Adding predicted values to the dataset
t= cor(Natural_MFE$problm_3, Natural_MFE$Fold)
ggplot(Natural_MFE, aes(Fold,problm_3)) + geom_point() + 
  ggtitle(label=(paste("R =", format(t, digits = 2)))) + 
  xlab("Measured Relative RNA Expression (log2)") + 
  ylab("Predicted Relative RNA Expression (log2)") + 
  geom_smooth(method = lm) +
  theme_science()

#read TE results
TE_results= read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/PARNA_RNA.txt")
TE_results_natural <- TE_results
#read grouped barcode identities
Grouped_barcode_identities = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Grouped_Barcode_identities.txt")
#filter RNA results with only those contianing Natural Let-7 binding sites
TE_results_natural <- dplyr::filter(TE_results_natural, RE_Identity %in% Grouped_barcode_identities$Natural)
#rename '7777' as 'Synthetic'
TE_results_natural$RE_Identity <- gsub("7777", "Synthetic", TE_results_natural$RE_Identity)
#plot TE
ggplot(TE_results_natural, aes(x = reorder(RE_Identity, -FOLD), FOLD)) + 
  xlab("") + ylab("Relative TE (log2)") +  
  geom_boxplot(outlier.shape = 1, notch =TRUE) + theme_science() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")
ggsave("Natural_TE.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#plot TE with invisible outliers
natural_TE <- ggplot(TE_results_natural, aes(x = reorder(RE_Identity, -FOLD), FOLD)) +  
  xlab("") + ylab("Relative TE (log2)") + 
  geom_boxplot(outlier.colour = NA, notch =TRUE) + theme_science() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")
#rescale above plot
sts_natural_2 <- boxplot.stats(TE_results_natural$FOLD)$stats
p2 = natural_TE + coord_cartesian(ylim = c(min(sts_natural_2)*1.8, max(sts_natural_2)))
p2
ggsave("Natural_TE_NO.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#assign an identifier for RNA or TE
TE_results_natural$Level <- "TE"
RNA_results_natural$Level <- "RNA"

#merge RNA and TE results
results_natural <- rbind(TE_results_natural, RNA_results_natural)
#plot fold change of RNA and TE
ggplot(results_natural, aes(reorder(RE_Identity, -FOLD), FOLD, fill=Level)) + 
  scale_fill_manual(values = c("gray", "white")) + 
  xlab("RE_Identity") + ylab("Fold Change (log2)") +
  geom_boxplot(outlier.shape = 1, notch =TRUE) + theme_science() + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank())
ggsave("natural_TE_RNA.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#plot fold change of RNA and TE with invisible outliers
natural_TE_RNA <- ggplot(results_natural, aes(reorder(RE_Identity, -FOLD), FOLD, fill=Level)) + 
  scale_fill_manual(values = c("gray", "white")) +
  xlab("RE Identity") + ylab("Fold Change (log2)") +
  geom_boxplot(outlier.colour = NA, notch =TRUE) + theme_science() + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), legend.title = element_text(size=14), legend.text = element_text(size=12), axis.text = element_text(size=12))

#rescale above plot
sts_5 <- boxplot.stats(results_natural$FOLD)$stats
p5 = natural_TE_RNA + coord_cartesian(ylim = c(min(sts_5),max(sts_5)*0.7))
p5

ggsave("natural_TE_RNA_NO.tiff", plot=last_plot(), width=5, height=5, units="in", dpi=300)
