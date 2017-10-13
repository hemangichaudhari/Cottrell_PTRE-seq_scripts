#open libraries
library(dplyr)
library(ggplot2)
library(extrafont)
library(gridExtra)
library(RGraphics)
font_import()

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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



#open sequencing results
RNA1_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L001_01_counts")
RNA2_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L001_02_counts")

RNA1_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L002_01_counts")
RNA2_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L002_02_counts")

RNA1_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L003_01_counts")
RNA2_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L003_02_counts")

RNA1_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L004_01_counts")
RNA2_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L004_02_counts")

#bind total RNA sequecing results for replicate 1
RNA1 = rbind(RNA1_L1, RNA1_L2, RNA1_L3, RNA1_L4)
RNA1 <- dplyr::group_by(RNA1, V2)
RNA1 <- dplyr::summarise(RNA1, Counts=sum(V1))
#bind total RNA sequecing results for replicate 2
RNA2 = rbind(RNA2_L1, RNA2_L2, RNA2_L3, RNA2_L4)
RNA2 <- dplyr::group_by(RNA2, V2)
RNA2 <- dplyr::summarise(RNA2, Counts=sum(V1))
#open sequencing results for 40s fraction
small_1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/Results_01_counts")
small_2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/Results_02_counts")



#change column names
colnames(RNA1) = c("Barcode", "Counts_RNA")
colnames(RNA2) = c("Barcode", "Counts_RNA")
colnames(small_1) = c("Counts_small", "Barcode")
colnames(small_2) = c("Counts_small", "Barcode")



#open "Barcode_Identities"
barcodes = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/Barcode_identities.txt")

#join "Barcode_Identities" with Sequencing Results
#merge sequencing results with barcodes
RNA1_id = merge(RNA1,barcodes, by ="Barcode")
RNA2_id = merge(RNA2,barcodes, by ="Barcode")

small_1_id = merge(small_1,barcodes, by ="Barcode")
small_2_id = merge(small_2,barcodes, by ="Barcode")




#Omit barcodes with no reads (NA)
RNA1_id = na.omit(RNA1_id)
RNA2_id = na.omit(RNA2_id)

small_1_id = na.omit(small_1_id)
small_2_id = na.omit(small_2_id)

#identify replicate in new column
small_1_id_edit <- data.frame(RE_Identity = small_1_id$RE_Identity, Barcode = small_1_id$Barcode, counts_small = small_1_id$Counts_small)
small_1_id_edit$Replicate <- 1
small_2_id_edit <- data.frame(RE_Identity = small_2_id$RE_Identity, Barcode = small_2_id$Barcode, counts_small = small_2_id$Counts_small)
small_2_id_edit$Replicate <- 2

#bind 40s data from each replicate
small_all <- rbind(small_1_id_edit, small_2_id_edit)
#rename some reporters
small_all$RE_Identity <- gsub("h", "A", small_all$RE_Identity)
small_all$RE_Identity <- gsub("7", "L", small_all$RE_Identity)
small_all$RE_Identity <- gsub("Lmer", "7mer", small_all$RE_Identity)
#write table for 40s sequencing
write.table(small_all, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/40SRNA_reads.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)




#Group barcodes by RE_Identity
RNA1_id<- dplyr::group_by(RNA1_id, RE_Identity)
RNA2_id <- dplyr::group_by(RNA2_id, RE_Identity)

small_1_id<- dplyr::group_by(small_1_id, RE_Identity)
small_2_id <- dplyr::group_by(small_1_id, RE_Identity)



#join cDNA and small by barcode then omit "na"

RNA1_small_1 = merge(RNA1_id, small_1_id, by ="Barcode")
RNA1_small_1_na = na.omit(RNA1_small_1)
RNA2_small_2 = merge(RNA2_id, small_2_id, by ="Barcode")
RNA2_small_2_na = na.omit(RNA2_small_2)


#normalize small_counts to cdNDA_counts
RNA1_small_1_na$normalized = RNA1_small_1_na$Counts_small /RNA1_small_1_na$Counts_RNA
RNA2_small_2_na$normalized = RNA2_small_2_na$Counts_small /RNA2_small_2_na$Counts_RNA


#make new data.frame with just barcode identities, and normalized values
RNA1_small_1_normalized = data.frame(RNA1_small_1_na$RE_Identity.x,RNA1_small_1_na$normalized)
colnames(RNA1_small_1_normalized) = c("RE_Identity","Normalized")
RNA2_small_2_normalized = data.frame(RNA2_small_2_na$RE_Identity.x,RNA2_small_2_na$normalized)
colnames(RNA2_small_2_normalized) = c("RE_Identity","Normalized")



#Define 'x' as the value for the control "BBBB"
x1 = subset(RNA1_small_1_normalized, RE_Identity=="BBBB")
x2 = subset(RNA2_small_2_normalized, RE_Identity=="BBBB")


#Add a new collumn with relative data for all RE and name it "Fold"
RNA1_small_1_relative <- dplyr::mutate(RNA1_small_1_normalized, log2(Normalized/median(x1$Normalized)))
colnames(RNA1_small_1_relative)[colnames(RNA1_small_1_relative)=="log2(Normalized/median(x1$Normalized))"] <- "FOLD"
RNA1_small_1__control <- RNA1_small_1_relative
RNA1_small_1__control$RE_Identity <- gsub("BBBB", "Control", RNA1_small_1__control$RE_Identity)

RNA2_small_2_relative <- dplyr::mutate(RNA2_small_2_normalized, log2(Normalized/median(x2$Normalized)))
colnames(RNA2_small_2_relative)[colnames(RNA2_small_2_relative)=="log2(Normalized/median(x2$Normalized))"] <- "FOLD"
RNA2_small_2__control <- RNA2_small_2_relative
RNA2_small_2__control$RE_Identity <- gsub("BBBB", "Control", RNA2_small_2__control$RE_Identity)



#group by RE_Identity and the summarise
RNA1_small_1_relative_grouped <- dplyr::group_by(RNA1_small_1_relative, RE_Identity) 
RNA1_small_1_relative_sum <- dplyr::summarise(RNA1_small_1_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1))

RNA2_small_2_relative_grouped <- dplyr::group_by(RNA2_small_2_relative, RE_Identity) 
RNA2_small_2_relative_sum <- dplyr::summarise(RNA2_small_2_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1))

#plot replicate to replicate comparisons
replicates_1_2 <- merge(RNA1_small_1_relative_sum, RNA2_small_2_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p1 <- ggplot(replicates_1_2, aes(x=Fold.x, y=Fold.y)) + xlab("Replicate 1") + ylab("Replicate 2") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  theme_science() + annotate(geom="text", x=1, y=-0.5, fontface=2, size=6, label="r = 0.926")
p1  
ggsave("small_replicates.tiff", height=4, width=4, units="in")
cor.test(replicates_1_2$Fold.x, replicates_1_2$Fold.y)


#bind replicate results
small_RNA <- rbind(RNA1_small_1_relative, RNA2_small_2_relative)
#Substitute the name for BBBB with Control
small_RNA$RE_Identity <- gsub("BBBB", "Control", small_RNA$RE_Identity)
#make a new column called "condensed" with all "B" removed
small_RNA$condensed <- gsub("B", "", small_RNA$RE_Identity)
#group by RE_Identity and the summarise
small_RNA_grouped <- dplyr::group_by(small_RNA, RE_Identity) 
small_RNA_relative_sum <- dplyr::summarise(small_RNA_grouped, Fold=median(FOLD), mad=mad(FOLD, constant=1), Fold_mean=mean(FOLD), sd=sd(FOLD))
small_RNA_relative_sum_level <- small_RNA_relative_sum
#assign group 40s
small_RNA_relative_sum_level$Group <- "40s"

#read TE and RNA results
TE_RNA_bound <- read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/TE_RNA_bound.txt")
#bind 40s and TE / RNA results
TE_RNA_40s <- rbind(small_RNA_relative_sum_level, TE_RNA_bound)

#define color palette
cbPalette <- c("#009E73", "#E69F00", "#56B4E9")
#plot TE, RNA and 40s
ggplot(TE_RNA_40s , aes(x=reorder(RE_Identity, -Fold), y=Fold, colour=Group)) + 
  scale_colour_manual(values= cbPalette, name="Level") +
  geom_point(alpha=0.5) + xlab("RE Identity") + 
  ylab("Fold Change (log2)") + theme_science() + 
  geom_hline(yintercept = 0, linetype="dashed") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("Boom_all.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#read hand-annotated TE and RNA table
annotated <- read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/TE_RNA_annotated.txt")
#read summarised TE and RNA data
TE_RNA_sum <- read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/TE_RNA.txt")
#merge annotations to TE and RNA data
TE_RNA_annotated <- merge(TE_RNA_sum, annotated, by="RE_Identity")
#rename a few regulatory elements
TE_RNA_annotated$group.y <- gsub("HuR", "ARE", TE_RNA_annotated$group.y)
TE_RNA_annotated$group.y <- gsub("Pumilio", "PRE", TE_RNA_annotated$group.y)
TE_RNA_annotated$group.y <- gsub("Smaug", "SRE", TE_RNA_annotated$group.y)
#assign legend order
legend_order <- c("Control", "Combination", "Let-7", "PRE", "ARE", "SRE")
#merge TE and RNA with 40s data
TE_RNA_40s_annotated <- merge(TE_RNA_annotated, small_RNA_relative_sum)
#define color palette
cbPalette <- c("#D55E00", "gray", "black", "#56B4E9", "#009E73", "#E69F00")
#plot TE vs 40s
ggplot(TE_RNA_40s_annotated, aes(x = Fold, y= Fold_TE, colour=group.y)) + 
  scale_colour_manual(name = "CRE", values=cbPalette, breaks = legend_order) + 
  xlab("Relative 40s Assoc. (log2)") + ylab("Relative TE (log2)") + coord_flip() +
  geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
  geom_point(alpha=0.5) + theme_science()
ggsave("TE_40s.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)
#plot RNA vs 40s
ggplot(TE_RNA_40s_annotated, aes(x = Fold, y= Fold_RNA, colour=group.y)) + 
  scale_colour_manual(name = "CRE", values=cbPalette, breaks = legend_order) + 
  xlab("Relative 40s Assoc. (log2)") + ylab("Relative Expression (log2)") + coord_flip() +
  geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
  geom_point(alpha=0.5) + theme_science()
ggsave("RNA_40s.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#group by Condensed and the summarise
small_RNA_grouped <- dplyr::group_by(small_RNA, condensed) 
small_RNA_condensed_sum <- dplyr::summarise(small_RNA_grouped, Fold_TE=median(FOLD), mad_TE=mad(FOLD, constant=1), Fold_TE_mean=mean(FOLD), sd_TE=sd(FOLD))

#write 40s data table
write.table(small_RNA, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/small_RNA.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)


