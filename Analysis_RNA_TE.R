#open libraries
library(dplyr)
library(ggplot2)
library(extrafont)
library(gridExtra)
library(RGraphics)
font_import()

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



#open sequencing results (PARNA =  polysome associated RNA)
RNA1_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L001_01_counts")
RNA2_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L001_02_counts")
RNA3_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L001_01_counts")
RNA4_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L001_02_counts")
PARNA1_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L001_03_counts")
PARNA2_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L001_04_counts")
PARNA3_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L001_03_counts")
PARNA4_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L001_04_counts")
Plasmid_1_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L001_05_counts")
Plasmid_2_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L001_05_counts")

RNA1_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L002_01_counts")
RNA2_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L002_02_counts")
RNA3_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L002_01_counts")
RNA4_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L002_02_counts")
PARNA1_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L002_03_counts")
PARNA2_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L002_04_counts")
PARNA3_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L002_03_counts")
PARNA4_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L002_04_counts")
Plasmid_1_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L002_05_counts")
Plasmid_2_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L002_05_counts")

RNA1_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L003_01_counts")
RNA2_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L003_02_counts")
RNA3_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L003_01_counts")
RNA4_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L003_02_counts")
PARNA1_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L003_03_counts")
PARNA2_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L003_04_counts")
PARNA3_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L003_03_counts")
PARNA4_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L003_04_counts")
Plasmid_1_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L003_05_counts")
Plasmid_2_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L003_05_counts")


RNA1_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L004_01_counts")
RNA2_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L004_02_counts")
RNA3_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L004_01_counts")
RNA4_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L004_02_counts")
PARNA1_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L004_03_counts")
PARNA2_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L004_04_counts")
PARNA3_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L004_03_counts")
PARNA4_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L004_04_counts")
Plasmid_1_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L004_05_counts")
Plasmid_2_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L004_05_counts")

#bind sequencing from each lane and find sum of counts for each barcode
RNA1 = rbind(RNA1_L1, RNA1_L2, RNA1_L3, RNA1_L4)
RNA1 <- dplyr::group_by(RNA1, V2)
RNA1 <- dplyr::summarise(RNA1, Counts=sum(V1))

RNA2 = rbind(RNA2_L1, RNA2_L2, RNA2_L3, RNA2_L4)
RNA2 <- dplyr::group_by(RNA2, V2)
RNA2 <- dplyr::summarise(RNA2, Counts=sum(V1))

RNA3 = rbind(RNA3_L1, RNA3_L2, RNA3_L3, RNA3_L4)
RNA3 <- dplyr::group_by(RNA3, V2)
RNA3 <- dplyr::summarise(RNA3, Counts=sum(V1))

RNA4 = rbind(RNA4_L1, RNA4_L2, RNA4_L3, RNA4_L4)
RNA4 <- dplyr::group_by(RNA4, V2)
RNA4 <- dplyr::summarise(RNA4, Counts=sum(V1))

PARNA1 = rbind(PARNA1_L1, PARNA1_L2, PARNA1_L3, PARNA1_L4)
PARNA1 <- dplyr::group_by(PARNA1, V2)
PARNA1 <- dplyr::summarise(PARNA1, Counts=sum(V1))

PARNA2 = rbind(PARNA2_L1, PARNA2_L2, PARNA2_L3, PARNA2_L4)
PARNA2 <- dplyr::group_by(PARNA2, V2)
PARNA2 <- dplyr::summarise(PARNA2, Counts=sum(V1))

PARNA3 = rbind(PARNA3_L1, PARNA3_L2, PARNA3_L3, PARNA3_L4)
PARNA3 <- dplyr::group_by(PARNA3, V2)
PARNA3 <- dplyr::summarise(PARNA3, Counts=sum(V1))

PARNA4 = rbind(PARNA4_L1, PARNA4_L2, PARNA4_L3, PARNA4_L4)
PARNA4 <- dplyr::group_by(PARNA4, V2)
PARNA4 <- dplyr::summarise(PARNA4, Counts=sum(V1))

Plasmid_1 = rbind(Plasmid_1_L1, Plasmid_1_L2, Plasmid_1_L3, Plasmid_1_L4)
Plasmid_1 <- dplyr::group_by(Plasmid_1, V2)
Plasmid_1 <- dplyr::summarise(Plasmid_1, Counts=sum(V1))

Plasmid_2 = rbind(Plasmid_2_L1, Plasmid_2_L2, Plasmid_2_L3, Plasmid_2_L4)
Plasmid_2 <- dplyr::group_by(Plasmid_2, V2)
Plasmid_2 <- dplyr::summarise(Plasmid_2, Counts=sum(V1))

#change column names
colnames(Plasmid_1) = c("Barcode", "Counts_Plasmid")
colnames(Plasmid_2) = c("Barcode", "Counts_Plasmid")
colnames(RNA1) = c("Barcode", "Counts_RNA")
colnames(RNA2) = c("Barcode", "Counts_RNA")
colnames(RNA3) = c("Barcode", "Counts_RNA")
colnames(RNA4) = c("Barcode", "Counts_RNA")
colnames(PARNA1) = c("Barcode", "Counts_PARNA")
colnames(PARNA2) = c("Barcode", "Counts_PARNA")
colnames(PARNA3) = c("Barcode", "Counts_PARNA")
colnames(PARNA4) = c("Barcode", "Counts_PARNA")


#open "Barcode_Identities"
barcodes = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Barcode_identities.txt")

#join "Barcode_Identities" with Sequencing Results
Plasmid_1_id = merge(Plasmid_1, barcodes, by="Barcode")
Plasmid_2_id = merge(Plasmid_2, barcodes, by="Barcode")

RNA1_id =merge(RNA1,barcodes, by ="Barcode")
RNA2_id =merge(RNA2,barcodes, by ="Barcode")
RNA3_id =merge(RNA3,barcodes, by ="Barcode")
RNA4_id =merge(RNA4,barcodes, by ="Barcode")

PARNA1_id =merge(PARNA1,barcodes, by ="Barcode")
PARNA2_id =merge(PARNA2,barcodes, by ="Barcode")
PARNA3_id =merge(PARNA3,barcodes, by ="Barcode")
PARNA4_id =merge(PARNA4,barcodes, by ="Barcode")

#Omit barcodes with no reads (NA)
Plasmid_1_id = na.omit(Plasmid_1_id)
Plasmid_2_id = na.omit(Plasmid_2_id)

#make new data.frame for plasmid sequencing and identify replicate in new column "Replicate"
Plasmid_1_id_edit <- data.frame(RE_Identity = Plasmid_1_id$RE_Identity, Barcode = Plasmid_1_id$Barcode, counts_Plasmid = Plasmid_1_id$Counts_Plasmid)
Plasmid_1_id_edit$Replicate <- 1
Plasmid_2_id_edit <- data.frame(RE_Identity = Plasmid_2_id$RE_Identity, Barcode = Plasmid_2_id$Barcode, counts_Plasmid = Plasmid_2_id$Counts_Plasmid)
Plasmid_2_id_edit$Replicate <- 2
#bind plasmid sequencing data
Plasmid_all <- rbind(Plasmid_1_id_edit, Plasmid_2_id_edit)
#rename h -> A (ARE) and 7 -> L (Let-7)
Plasmid_all$RE_Identity <- gsub("h", "A", Plasmid_all$RE_Identity)
Plasmid_all$RE_Identity <- gsub("7", "L", Plasmid_all$RE_Identity)
Plasmid_all$RE_Identity <- gsub("Lmer", "7mer", Plasmid_all$RE_Identity)
#write.table of plasmid sequencing barcode counts
write.table(Plasmid_all, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Plasmid_reads.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

#Omit barcodes with no reads (NA)
RNA1_id = na.omit(RNA1_id)
RNA2_id = na.omit(RNA2_id)
RNA3_id = na.omit(RNA3_id)
RNA4_id = na.omit(RNA4_id)

#make new data.frame for sequencing and identify replicate in new column "Replicate"
RNA1_id_edit <- data.frame(RE_Identity = RNA1_id$RE_Identity, Barcode = RNA1_id$Barcode, counts_RNA = RNA1_id$Counts_RNA)
RNA1_id_edit$Replicate <- 1
RNA2_id_edit <- data.frame(RE_Identity = RNA2_id$RE_Identity, Barcode = RNA2_id$Barcode, counts_RNA = RNA2_id$Counts_RNA)
RNA2_id_edit$Replicate <- 2
RNA3_id_edit <- data.frame(RE_Identity = RNA3_id$RE_Identity, Barcode = RNA3_id$Barcode, counts_RNA = RNA3_id$Counts_RNA)
RNA3_id_edit$Replicate <- 3
RNA4_id_edit <- data.frame(RE_Identity = RNA4_id$RE_Identity, Barcode = RNA4_id$Barcode, counts_RNA = RNA4_id$Counts_RNA)
RNA4_id_edit$Replicate <- 4
#bind RNA sequenicng results
RNA_all <- rbind(RNA1_id_edit, RNA2_id_edit, RNA3_id_edit, RNA4_id_edit)
#rename h -> A (ARE) and 7 -> L (Let-7)
RNA_all$RE_Identity <- gsub("h", "A", RNA_all$RE_Identity)
RNA_all$RE_Identity <- gsub("7", "L", RNA_all$RE_Identity)
RNA_all$RE_Identity <- gsub("Lmer", "7mer", RNA_all$RE_Identity)
#write.table of RNA sequencing barcode counts
write.table(RNA_all, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/RNA_reads.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

#Omit barcodes with no reads (NA)
PARNA1_id = na.omit(PARNA1_id)
PARNA2_id = na.omit(PARNA2_id)
PARNA3_id = na.omit(PARNA3_id)
PARNA4_id = na.omit(PARNA4_id)

#make new data.frame for sequencing and identify replicate in new column "Replicate"
PARNA1_id_edit <- data.frame(RE_Identity = PARNA1_id$RE_Identity, Barcode = PARNA1_id$Barcode, counts_PARNA = PARNA1_id$Counts_PARNA)
PARNA1_id_edit$Replicate <- 1
PARNA2_id_edit <- data.frame(RE_Identity = PARNA2_id$RE_Identity, Barcode = PARNA2_id$Barcode, counts_PARNA = PARNA2_id$Counts_PARNA)
PARNA2_id_edit$Replicate <- 2
PARNA3_id_edit <- data.frame(RE_Identity = PARNA3_id$RE_Identity, Barcode = PARNA3_id$Barcode, counts_PARNA = PARNA3_id$Counts_PARNA)
PARNA3_id_edit$Replicate <- 3
PARNA4_id_edit <- data.frame(RE_Identity = PARNA4_id$RE_Identity, Barcode = PARNA4_id$Barcode, counts_PARNA = PARNA4_id$Counts_PARNA)
PARNA4_id_edit$Replicate <- 4
#bind PARNA sequenicng results
PARNA_all <- rbind(PARNA1_id_edit, PARNA2_id_edit, PARNA3_id_edit, PARNA4_id_edit)
#rename h -> A (ARE) and 7 -> L (Let-7)
PARNA_all$RE_Identity <- gsub("h", "A", PARNA_all$RE_Identity)
PARNA_all$RE_Identity <- gsub("7", "L", PARNA_all$RE_Identity)
PARNA_all$RE_Identity <- gsub("Lmer", "7mer", PARNA_all$RE_Identity)
#write.table of PARNA sequencing barcode counts
write.table(PARNA_all, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/PARNA_reads.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)


#Omit barcodes with less than 10 reads

Plasmid_1_id = subset(Plasmid_1_id, Counts_Plasmid > 10)
Plasmid_2_id = subset(Plasmid_2_id, Counts_Plasmid > 10)

RNA1_id = subset(RNA1_id, Counts_RNA >10)
RNA2_id = subset(RNA2_id, Counts_RNA >10)
RNA3_id = subset(RNA3_id, Counts_RNA >10)
RNA4_id = subset(RNA4_id, Counts_RNA >10)

PARNA1_id = subset(PARNA1_id, Counts_PARNA >10)
PARNA2_id = subset(PARNA2_id, Counts_PARNA >10)
PARNA3_id = subset(PARNA3_id, Counts_PARNA >10)
PARNA4_id = subset(PARNA4_id, Counts_PARNA >10)

#merge plasmid sequencing results
Plasmid <- merge(Plasmid_1_id, Plasmid_2_id, by.x = "Barcode", by.y = "Barcode")
#find average barcode counts for each sequencing run for plasmid
Plasmid$Counts_avg <- (Plasmid$Counts_Plasmid.x + Plasmid$Counts_Plasmid.y)/2

#Group barcodes by RE_Identity
Plasmid  <- dplyr::group_by(Plasmid , RE_Identity.x) 

RNA1_id<- dplyr::group_by(RNA1_id, RE_Identity)
RNA2_id <- dplyr::group_by(RNA2_id, RE_Identity)
RNA3_id <- dplyr::group_by(RNA3_id, RE_Identity)
RNA4_id <- dplyr::group_by(RNA4_id, RE_Identity)

PARNA1_id<- dplyr::group_by(PARNA1_id, RE_Identity)
PARNA2_id <- dplyr::group_by(PARNA2_id, RE_Identity)
PARNA3_id <- dplyr::group_by(PARNA3_id, RE_Identity)
PARNA4_id <- dplyr::group_by(PARNA4_id, RE_Identity)

#plot hist of barcode counts
hist(log10(Plasmid$Counts_avg))

hist(log10(RNA1_id$Counts_RNA))
hist(log10(RNA2_id$Counts_RNA))
hist(log10(RNA3_id$Counts_RNA))
hist(log10(RNA4_id$Counts_RNA))

hist(log10(PARNA1_id$Counts_PARNA))
hist(log10(PARNA2_id$Counts_PARNA))
hist(log10(PARNA3_id$Counts_PARNA))
hist(log10(PARNA4_id$Counts_PARNA))

#join cDNA and plasmid by barcode then omit "na"

RNA1_plasmid = merge(RNA1_id, Plasmid, by ="Barcode")
RNA1_plasmid_na = na.omit(RNA1_plasmid)
RNA2_plasmid = merge(RNA2_id,Plasmid, by ="Barcode")
RNA2_plasmid_na = na.omit(RNA2_plasmid)
RNA3_plasmid = merge(RNA3_id,Plasmid, by ="Barcode")
RNA3_plasmid_na = na.omit(RNA3_plasmid)
RNA4_plasmid = merge(RNA4_id,Plasmid, by ="Barcode")
RNA4_plasmid_na = na.omit(RNA4_plasmid)


#normalize cDNA_counts to plasmid_counts
RNA1_plasmid_na$normalized = RNA1_plasmid_na$Counts_RNA/RNA1_plasmid_na$Counts_avg
RNA2_plasmid_na$normalized = RNA2_plasmid_na$Counts_RNA/RNA2_plasmid_na$Counts_avg
RNA3_plasmid_na$normalized = RNA3_plasmid_na$Counts_RNA/RNA3_plasmid_na$Counts_avg
RNA4_plasmid_na$normalized = RNA4_plasmid_na$Counts_RNA/RNA4_plasmid_na$Counts_avg

#make new data.frame with just barcode identities, and normalized values
RNA1_plasmid_normalized = data.frame(RNA1_plasmid_na$RE_Identity.x, RNA1_plasmid_na$normalized)
colnames(RNA1_plasmid_normalized) = c("RE_Identity","Normalized")
RNA2_plasmid_normalized = data.frame(RNA2_plasmid_na$RE_Identity.x, RNA2_plasmid_na$normalized)
colnames(RNA2_plasmid_normalized) = c("RE_Identity","Normalized")
RNA3_plasmid_normalized = data.frame(RNA3_plasmid_na$RE_Identity.x, RNA3_plasmid_na$normalized)
colnames(RNA3_plasmid_normalized) = c("RE_Identity","Normalized")
RNA4_plasmid_normalized = data.frame(RNA4_plasmid_na$RE_Identity.x, RNA4_plasmid_na$normalized)
colnames(RNA4_plasmid_normalized) = c("RE_Identity","Normalized")

#Make subset of data.frame for control "BBBB"
x1 = subset(RNA1_plasmid_normalized, RE_Identity=="BBBB")
x2 = subset(RNA2_plasmid_normalized, RE_Identity=="BBBB")
x3 = subset(RNA3_plasmid_normalized, RE_Identity=="BBBB")
x4 = subset(RNA4_plasmid_normalized, RE_Identity=="BBBB")

#Find 'Fold' change in expression
RNA1_plasmid_relative <- dplyr::mutate(RNA1_plasmid_normalized, log2(Normalized/median(x1$Normalized)))
colnames(RNA1_plasmid_relative)[colnames(RNA1_plasmid_relative)=="log2(Normalized/median(x1$Normalized))"] <- "FOLD"
RNA1_plasmid_control <- RNA1_plasmid_relative
RNA1_plasmid_control$RE_Identity <- gsub("BBBB", "Control", RNA1_plasmid_control$RE_Identity)

RNA2_plasmid_relative <- dplyr::mutate(RNA2_plasmid_normalized, log2(Normalized/median(x2$Normalized)))
colnames(RNA2_plasmid_relative)[colnames(RNA2_plasmid_relative)=="log2(Normalized/median(x2$Normalized))"] <- "FOLD"
RNA2_plasmid_control <- RNA2_plasmid_relative
RNA2_plasmid_control$RE_Identity <- gsub("BBBB", "Control", RNA2_plasmid_control$RE_Identity)

RNA3_plasmid_relative <- dplyr::mutate(RNA3_plasmid_normalized, log2(Normalized/median(x3$Normalized)))
colnames(RNA3_plasmid_relative)[colnames(RNA3_plasmid_relative)=="log2(Normalized/median(x3$Normalized))"] <- "FOLD"
RNA3_plasmid_control <- RNA3_plasmid_relative
RNA3_plasmid_control$RE_Identity <- gsub("BBBB", "Control", RNA3_plasmid_control$RE_Identity)

RNA4_plasmid_relative <- dplyr::mutate(RNA4_plasmid_normalized, log2(Normalized/median(x4$Normalized)))
colnames(RNA4_plasmid_relative)[colnames(RNA4_plasmid_relative)=="log2(Normalized/median(x4$Normalized))"] <- "FOLD"
RNA4_plasmid_control <- RNA4_plasmid_relative
RNA4_plasmid_control$RE_Identity <- gsub("BBBB", "Control", RNA4_plasmid_control$RE_Identity)


#group by RE_Identity and the summarise

RNA1_plasmid_relative_grouped <- dplyr::group_by(RNA1_plasmid_relative, RE_Identity) 
RNA1_plasmid_relative_sum <- dplyr::summarise(RNA1_plasmid_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), stdev=sd(FOLD), se=se(FOLD), count=length(FOLD))
RNA1_plasmid_relative_sum$Replicate <- "Replicate_1"

RNA2_plasmid_relative_grouped <- dplyr::group_by(RNA2_plasmid_relative, RE_Identity) 
RNA2_plasmid_relative_sum <- dplyr::summarise(RNA2_plasmid_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), stdev=sd(FOLD), se=se(FOLD), count=length(FOLD))
RNA2_plasmid_relative_sum$Replicate <- "Replicate_2"

RNA3_plasmid_relative_grouped <- dplyr::group_by(RNA3_plasmid_relative, RE_Identity) 
RNA3_plasmid_relative_sum <- dplyr::summarise(RNA3_plasmid_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), stdev=sd(FOLD), se=se(FOLD), count=length(FOLD))
RNA3_plasmid_relative_sum$Replicate <- "Replicate_3"

RNA4_plasmid_relative_grouped <- dplyr::group_by(RNA4_plasmid_relative, RE_Identity) 
RNA4_plasmid_relative_sum <- dplyr::summarise(RNA4_plasmid_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), stdev=sd(FOLD), se=se(FOLD), count=length(FOLD))
RNA4_plasmid_relative_sum$Replicate <- "Replicate_4"


#plot replicate to replicate comparisons of 'Fold Change'
replicates_1_2 <- merge(RNA1_plasmid_relative_sum , RNA2_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p1 <- ggplot(replicates_1_2, aes(x=Fold.x, y=Fold.y)) + xlab("Replicate 1") + ylab("Replicate 2") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p1  
cor.test(replicates_1_2$Fold.x, replicates_1_2$Fold.y)

replicates_1_3 <- merge(RNA1_plasmid_relative_sum , RNA3_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p2 <- ggplot(replicates_1_3, aes(x=Fold.x, y=Fold.y)) + xlab("") + ylab("Replicate 3") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p2 
cor.test(replicates_1_3$Fold.x, replicates_1_3$Fold.y)

replicates_1_4 <- merge(RNA1_plasmid_relative_sum , RNA4_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p3 <- ggplot(replicates_1_4, aes(x=Fold.x, y=Fold.y)) + xlab("") + ylab("Replicate 4") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p3
cor.test(replicates_1_4$Fold.x, replicates_1_4$Fold.y)

replicates_2_3 <- merge(RNA2_plasmid_relative_sum , RNA3_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p4 <- ggplot(replicates_2_3, aes(x=Fold.x, y=Fold.y)) + xlab("Replicate 2") + ylab("") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p4 
cor.test(replicates_2_3$Fold.x, replicates_2_3$Fold.y)

replicates_2_4 <- merge(RNA2_plasmid_relative_sum , RNA4_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p5 <- ggplot(replicates_2_4, aes(x=Fold.x, y=Fold.y)) + xlab("") + ylab("") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p5 
cor.test(replicates_2_4$Fold.x, replicates_2_4$Fold.y)

replicates_3_4 <- merge(RNA3_plasmid_relative_sum , RNA4_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p6 <- ggplot(replicates_3_4, aes(x=Fold.x, y=Fold.y)) + xlab("Replicate 3") + ylab("") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p6
cor.test(replicates_3_4$Fold.x, replicates_3_4$Fold.y)

pb <- textGrob("")

#determine pearson correlation coeff. for replicate to replicate comparisons
m_rna <- rbind(cor(replicates_1_2$Fold.x, replicates_1_2$Fold.y), 
               cor(replicates_1_3$Fold.x, replicates_1_3$Fold.y),
               cor(replicates_1_4$Fold.x, replicates_1_4$Fold.y),
               cor(replicates_2_3$Fold.x, replicates_2_3$Fold.y),
               cor(replicates_2_4$Fold.x, replicates_2_4$Fold.y),
               cor(replicates_3_4$Fold.x, replicates_3_4$Fold.y))
#labels for comparisons
Comparison <- c("Replicate 1 x 2","Replicate 1 x 3", "Replicate 1 x 4", "Replicate 2 x 3", "Replicate 2 x 4", "Replicate 3 x 4")
#data frame for pearson correlation coeff. for replicate to replicate comparisons
mean_rna <- data.frame(Comparison, m_rna)
colnames(mean_rna) <- c("Comparison", "Pearson Correlation")
#make table graphic for replicate to replicate comparisons
px <- tableGrob(mean_rna)
#arrange tables and graphs
grid.arrange(p3, p5, p6, p2, p4, pb, p1, pb, px, ncol=3)
g_rna <- arrangeGrob(p3, p5, p6, p2, p4, pb, p1, pb, px, ncol=3)
ggsave("replicates_rna.tiff", g_rna, height=8, width=10, dpi=300, units="in")

#join cDNA and plasmid by barcode then omit "na"

PARNA1_plasmid = merge(PARNA1_id, Plasmid, by ="Barcode")
PARNA1_plasmid_na = na.omit(PARNA1_plasmid)
PARNA2_plasmid = merge(PARNA2_id,Plasmid, by ="Barcode")
PARNA2_plasmid_na = na.omit(PARNA2_plasmid)
PARNA3_plasmid = merge(PARNA3_id,Plasmid, by ="Barcode")
PARNA3_plasmid_na = na.omit(PARNA3_plasmid)
PARNA4_plasmid = merge(PARNA4_id,Plasmid, by ="Barcode")
PARNA4_plasmid_na = na.omit(PARNA4_plasmid)

#normalize cDNA_counts to plasmid_counts
PARNA1_plasmid_na$normalized = PARNA1_plasmid_na$Counts_PARNA/PARNA1_plasmid_na$Counts_avg
PARNA2_plasmid_na$normalized = PARNA2_plasmid_na$Counts_PARNA/PARNA2_plasmid_na$Counts_avg
PARNA3_plasmid_na$normalized = PARNA3_plasmid_na$Counts_PARNA/PARNA3_plasmid_na$Counts_avg
PARNA4_plasmid_na$normalized = PARNA4_plasmid_na$Counts_PARNA/PARNA4_plasmid_na$Counts_avg

#make new data.frame with just barcode identities, and normalized values
PARNA1_plasmid_normalized = data.frame(PARNA1_plasmid_na$RE_Identity.x, PARNA1_plasmid_na$normalized)
colnames(PARNA1_plasmid_normalized) = c("RE_Identity","Normalized")
PARNA2_plasmid_normalized = data.frame(PARNA2_plasmid_na$RE_Identity.x, PARNA2_plasmid_na$normalized)
colnames(PARNA2_plasmid_normalized) = c("RE_Identity","Normalized")
PARNA3_plasmid_normalized = data.frame(PARNA3_plasmid_na$RE_Identity.x, PARNA3_plasmid_na$normalized)
colnames(PARNA3_plasmid_normalized) = c("RE_Identity","Normalized")
PARNA4_plasmid_normalized = data.frame(PARNA4_plasmid_na$RE_Identity.x, PARNA4_plasmid_na$normalized)
colnames(PARNA4_plasmid_normalized) = c("RE_Identity","Normalized")


#Define 'x' as the value for the control "BBBB"
y1 = subset(PARNA1_plasmid_normalized, RE_Identity=="BBBB")
y2 = subset(PARNA2_plasmid_normalized, RE_Identity=="BBBB")
y3 = subset(PARNA3_plasmid_normalized, RE_Identity=="BBBB")
y4 = subset(PARNA4_plasmid_normalized, RE_Identity=="BBBB")

#Find 'Fold' change of polysome normalized by plasmid
PARNA1_plasmid_relative <- dplyr::mutate(PARNA1_plasmid_normalized, log2(Normalized/median(y1$Normalized)))
colnames(PARNA1_plasmid_relative)[colnames(PARNA1_plasmid_relative)=="log2(Normalized/median(y1$Normalized))"] <- "FOLD"
PARNA1_plasmid_control <- PARNA1_plasmid_relative
PARNA1_plasmid_control$RE_Identity <- gsub("BBBB", "Control", PARNA1_plasmid_control$RE_Identity)

PARNA2_plasmid_relative <- dplyr::mutate(PARNA2_plasmid_normalized, log2(Normalized/median(y2$Normalized)))
colnames(PARNA2_plasmid_relative)[colnames(PARNA2_plasmid_relative)=="log2(Normalized/median(y2$Normalized))"] <- "FOLD"
PARNA2_plasmid_control <- PARNA2_plasmid_relative
PARNA2_plasmid_control$RE_Identity <- gsub("BBBB", "Control", PARNA2_plasmid_control$RE_Identity)

PARNA3_plasmid_relative <- dplyr::mutate(PARNA3_plasmid_normalized, log2(Normalized/median(y3$Normalized)))
colnames(PARNA3_plasmid_relative)[colnames(PARNA3_plasmid_relative)=="log2(Normalized/median(y3$Normalized))"] <- "FOLD"
PARNA3_plasmid_control <- PARNA3_plasmid_relative
PARNA3_plasmid_control$RE_Identity <- gsub("BBBB", "Control", PARNA3_plasmid_control$RE_Identity)

PARNA4_plasmid_relative <- dplyr::mutate(PARNA4_plasmid_normalized, log2(Normalized/median(y4$Normalized)))
colnames(PARNA4_plasmid_relative)[colnames(PARNA4_plasmid_relative)=="log2(Normalized/median(y4$Normalized))"] <- "FOLD"
PARNA4_plasmid_control <- PARNA4_plasmid_relative
PARNA4_plasmid_control$RE_Identity <- gsub("BBBB", "Control", PARNA4_plasmid_control$RE_Identity)


#group by RE_Identity and the summarise
PARNA1_plasmid_relative_grouped <- dplyr::group_by(PARNA1_plasmid_relative, RE_Identity) 
PARNA1_plasmid_relative_sum <- dplyr::summarise(PARNA1_plasmid_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), count=length(FOLD))

PARNA2_plasmid_relative_grouped <- dplyr::group_by(PARNA2_plasmid_relative, RE_Identity) 
PARNA2_plasmid_relative_sum <- dplyr::summarise(PARNA2_plasmid_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), count=length(FOLD))

PARNA3_plasmid_relative_grouped <- dplyr::group_by(PARNA3_plasmid_relative, RE_Identity) 
PARNA3_plasmid_relative_sum <- dplyr::summarise(PARNA3_plasmid_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), count=length(FOLD))

PARNA4_plasmid_relative_grouped <- dplyr::group_by(PARNA4_plasmid_relative, RE_Identity) 
PARNA4_plasmid_relative_sum <- dplyr::summarise(PARNA4_plasmid_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), count=length(FOLD))



#plot replicate to replicate comparisons of 'Fold Change'
replicates_1_2 <- merge(PARNA1_plasmid_relative_sum , PARNA2_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p1 <- ggplot(replicates_1_2, aes(x=Fold.x, y=Fold.y)) + xlab("Replicate 1") + ylab("Replicate 2") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p1  
cor.test(replicates_1_2$Fold.x, replicates_1_2$Fold.y)

replicates_1_3 <- merge(PARNA1_plasmid_relative_sum , PARNA3_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p2 <- ggplot(replicates_1_3, aes(x=Fold.x, y=Fold.y)) + xlab("") + ylab("Replicate 3") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p2 
cor.test(replicates_1_3$Fold.x, replicates_1_3$Fold.y)

replicates_1_4 <- merge(PARNA1_plasmid_relative_sum , PARNA4_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p3 <- ggplot(replicates_1_4, aes(x=Fold.x, y=Fold.y)) + xlab("") + ylab("Replicate 4") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p3
cor.test(replicates_1_4$Fold.x, replicates_1_4$Fold.y)

replicates_2_3 <- merge(PARNA2_plasmid_relative_sum , PARNA3_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p4 <- ggplot(replicates_2_3, aes(x=Fold.x, y=Fold.y)) + xlab("Replicate 2") + ylab("") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p4 
cor.test(replicates_2_3$Fold.x, replicates_2_3$Fold.y)

replicates_2_4 <- merge(PARNA2_plasmid_relative_sum , PARNA4_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p5 <- ggplot(replicates_2_4, aes(x=Fold.x, y=Fold.y)) + xlab("") + ylab("") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p5 
cor.test(replicates_2_4$Fold.x, replicates_2_4$Fold.y)

replicates_3_4 <- merge(PARNA3_plasmid_relative_sum , PARNA4_plasmid_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p6 <- ggplot(replicates_3_4, aes(x=Fold.x, y=Fold.y)) + xlab("Replicate 3") + ylab("") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p6
cor.test(replicates_3_4$Fold.x, replicates_3_4$Fold.y)

pb <- textGrob("")
#determine pearson correlation coeff. for replicate to replicate comparisons
m_PP <- rbind(cor(replicates_1_2$Fold.x, replicates_1_2$Fold.y), 
              cor(replicates_1_3$Fold.x, replicates_1_3$Fold.y),
              cor(replicates_1_4$Fold.x, replicates_1_4$Fold.y),
              cor(replicates_2_3$Fold.x, replicates_2_3$Fold.y),
              cor(replicates_2_4$Fold.x, replicates_2_4$Fold.y),
              cor(replicates_3_4$Fold.x, replicates_3_4$Fold.y))

#labels for comparisons
Comparison <- c("Replicate 1 x 2","Replicate 1 x 3", "Replicate 1 x 4", "Replicate 2 x 3", "Replicate 2 x 4", "Replicate 3 x 4")
#data frame for pearson correlation coeff. for replicate to replicate comparisons
mean_PP <- data.frame(Comparison, m_PP)
colnames(mean_PP) <- c("Comparison", "Pearson Correlation")
#make table graphic for replicate to replicate comparisons
px <- tableGrob(mean_PP)
#arrange tables and graphs
grid.arrange(p3, p5, p6, p2, p4, pb, p1, pb, px, ncol=3)
g_PP <- arrangeGrob(p3, p5, p6, p2, p4, pb, p1, pb, px, ncol=3)
ggsave("replicates_PARNA.tiff", g_PP, height=8, width=10, dpi=300, units="in")

#join cDNA and plasmid by barcode then omit "na"

PARNA1_RNA = merge(PARNA1_id, RNA1_id, by ="Barcode")
PARNA1_RNA_na = na.omit(PARNA1_RNA)
PARNA2_RNA = merge(PARNA2_id,RNA2_id, by ="Barcode")
PARNA2_RNA_na = na.omit(PARNA2_RNA)
PARNA3_RNA = merge(PARNA3_id,RNA3_id, by ="Barcode")
PARNA3_RNA_na = na.omit(PARNA3_RNA)
PARNA4_RNA = merge(PARNA4_id,RNA4_id, by ="Barcode")
PARNA4_RNA_na = na.omit(PARNA4_RNA)

#normalize cDNA_counts to RNA_counts
PARNA1_RNA_na$normalized = PARNA1_RNA_na$Counts_PARNA/PARNA1_RNA_na$Counts_RNA
PARNA2_RNA_na$normalized = PARNA2_RNA_na$Counts_PARNA/PARNA2_RNA_na$Counts_RNA
PARNA3_RNA_na$normalized = PARNA3_RNA_na$Counts_PARNA/PARNA3_RNA_na$Counts_RNA
PARNA4_RNA_na$normalized = PARNA4_RNA_na$Counts_PARNA/PARNA4_RNA_na$Counts_RNA

ggplot(PARNA1_RNA_na, aes(Counts_RNA, normalized)) + geom_point()

#make new data.frame with just barcode identities, and normalized values
PARNA1_RNA_normalized = data.frame(PARNA1_RNA_na$RE_Identity.x, PARNA1_RNA_na$normalized)
colnames(PARNA1_RNA_normalized) = c("RE_Identity","Normalized")
PARNA2_RNA_normalized = data.frame(PARNA2_RNA_na$RE_Identity.x, PARNA2_RNA_na$normalized)
colnames(PARNA2_RNA_normalized) = c("RE_Identity","Normalized")
PARNA3_RNA_normalized = data.frame(PARNA3_RNA_na$RE_Identity.x, PARNA3_RNA_na$normalized)
colnames(PARNA3_RNA_normalized) = c("RE_Identity","Normalized")
PARNA4_RNA_normalized = data.frame(PARNA4_RNA_na$RE_Identity.x, PARNA4_RNA_na$normalized)
colnames(PARNA4_RNA_normalized) = c("RE_Identity","Normalized")


#Define 'x' as the value for the control "BBBB"
t1 = subset(PARNA1_RNA_normalized, RE_Identity=="BBBB")
t2 = subset(PARNA2_RNA_normalized, RE_Identity=="BBBB")
t3 = subset(PARNA3_RNA_normalized, RE_Identity=="BBBB")
t4 = subset(PARNA4_RNA_normalized, RE_Identity=="BBBB")

#Determine 'Fold' change of TE
PARNA1_RNA_relative <- dplyr::mutate(PARNA1_RNA_normalized, log2(Normalized/median(t1$Normalized)))
colnames(PARNA1_RNA_relative)[colnames(PARNA1_RNA_relative)=="log2(Normalized/median(t1$Normalized))"] <- "FOLD"
PARNA1_RNA_control <- PARNA1_RNA_relative
PARNA1_RNA_control$RE_Identity <- gsub("BBBB", "Control", PARNA1_RNA_control$RE_Identity)

PARNA2_RNA_relative <- dplyr::mutate(PARNA2_RNA_normalized, log2(Normalized/median(t2$Normalized)))
colnames(PARNA2_RNA_relative)[colnames(PARNA2_RNA_relative)=="log2(Normalized/median(t2$Normalized))"] <- "FOLD"
PARNA2_RNA_control <- PARNA2_RNA_relative
PARNA2_RNA_control$RE_Identity <- gsub("BBBB", "Control", PARNA2_RNA_control$RE_Identity)

PARNA3_RNA_relative <- dplyr::mutate(PARNA3_RNA_normalized, log2(Normalized/median(t3$Normalized)))
colnames(PARNA3_RNA_relative)[colnames(PARNA3_RNA_relative)=="log2(Normalized/median(t3$Normalized))"] <- "FOLD"
PARNA3_RNA_control <- PARNA3_RNA_relative
PARNA3_RNA_control$RE_Identity <- gsub("BBBB", "Control", PARNA3_RNA_control$RE_Identity)

PARNA4_RNA_relative <- dplyr::mutate(PARNA4_RNA_normalized, log2(Normalized/median(t4$Normalized)))
colnames(PARNA4_RNA_relative)[colnames(PARNA4_RNA_relative)=="log2(Normalized/median(t4$Normalized))"] <- "FOLD"
PARNA4_RNA_control <- PARNA2_RNA_relative
PARNA4_RNA_control$RE_Identity <- gsub("BBBB", "Control", PARNA2_RNA_control$RE_Identity)


#group by RE_Identity and the summarise
PARNA1_RNA_relative_grouped <- dplyr::group_by(PARNA1_RNA_relative, RE_Identity) 
PARNA1_RNA_relative_sum <- dplyr::summarise(PARNA1_RNA_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), count=length(FOLD))

PARNA2_RNA_relative_grouped <- dplyr::group_by(PARNA2_RNA_relative, RE_Identity) 
PARNA2_RNA_relative_sum <- dplyr::summarise(PARNA2_RNA_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), count=length(FOLD))

PARNA3_RNA_relative_grouped <- dplyr::group_by(PARNA3_RNA_relative, RE_Identity) 
PARNA3_RNA_relative_sum <- dplyr::summarise(PARNA3_RNA_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), count=length(FOLD))

PARNA4_RNA_relative_grouped <- dplyr::group_by(PARNA4_RNA_relative, RE_Identity) 
PARNA4_RNA_relative_sum <- dplyr::summarise(PARNA4_RNA_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1), count=length(FOLD))

#plot replicate to replicate comparisons
replicates_1_2 <- merge(PARNA1_RNA_relative_sum , PARNA2_RNA_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p1 <- ggplot(replicates_1_2, aes(x=Fold.x, y=Fold.y)) + xlab("Replicate 1") + ylab("Replicate 2") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p1  
cor.test(replicates_1_2$Fold.x, replicates_1_2$Fold.y)

replicates_1_3 <- merge(PARNA1_RNA_relative_sum , PARNA3_RNA_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p2 <- ggplot(replicates_1_3, aes(x=Fold.x, y=Fold.y)) + xlab("") + ylab("Replicate 3") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p2 
cor.test(replicates_1_3$Fold.x, replicates_1_3$Fold.y)

replicates_1_4 <- merge(PARNA1_RNA_relative_sum , PARNA4_RNA_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p3 <- ggplot(replicates_1_4, aes(x=Fold.x, y=Fold.y)) + xlab("") + ylab("Replicate 4") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p3
cor.test(replicates_1_4$Fold.x, replicates_1_4$Fold.y)

replicates_2_3 <- merge(PARNA2_RNA_relative_sum , PARNA3_RNA_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p4 <- ggplot(replicates_2_3, aes(x=Fold.x, y=Fold.y)) + xlab("Replicate 2") + ylab("") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p4 
cor.test(replicates_2_3$Fold.x, replicates_2_3$Fold.y)

replicates_2_4 <- merge(PARNA2_RNA_relative_sum , PARNA4_RNA_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p5 <- ggplot(replicates_2_4, aes(x=Fold.x, y=Fold.y)) + xlab("") + ylab("") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p5 
cor.test(replicates_2_4$Fold.x, replicates_2_4$Fold.y)

replicates_3_4 <- merge(PARNA3_RNA_relative_sum , PARNA4_RNA_relative_sum , by.x = "RE_Identity", by.y= "RE_Identity")
p6 <- ggplot(replicates_3_4, aes(x=Fold.x, y=Fold.y)) + xlab("Replicate 3") + ylab("") +
  geom_point()+ geom_abline() + stat_smooth(method="lm") + theme_science() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
p6
cor.test(replicates_3_4$Fold.x, replicates_3_4$Fold.y)

pb <- textGrob("")
#find pearson correlation coeffecient for replicate to replicate comparisons
m_te <- rbind(cor(replicates_1_2$Fold.x, replicates_1_2$Fold.y), 
              cor(replicates_1_3$Fold.x, replicates_1_3$Fold.y),
              cor(replicates_1_4$Fold.x, replicates_1_4$Fold.y),
              cor(replicates_2_3$Fold.x, replicates_2_3$Fold.y),
              cor(replicates_2_4$Fold.x, replicates_2_4$Fold.y),
              cor(replicates_3_4$Fold.x, replicates_3_4$Fold.y))
#labels for comparisons
Comparison <- c("Replicate 1 x 2","Replicate 1 x 3", "Replicate 1 x 4", "Replicate 2 x 3", "Replicate 2 x 4", "Replicate 3 x 4")
#make data.frame for pearsons correlation coefficients
mean_te <- data.frame(Comparison, m_te)
colnames(mean_te) <- c("Comparison", "Pearson Correlation")
#make table graphic for pearson coeff.
px <- tableGrob(mean_te)
#arrange plots
grid.arrange(p3, p5, p6, p2, p4, pb, p1, pb, px, ncol=3)
g_TE <- arrangeGrob(p3, p5, p6, p2, p4, pb, p1, pb, px, ncol=3)
ggsave("replicates_TE.tiff", g_TE, height=8, width=10, dpi=300, units="in")

#bind replicate results for TE
PARNA_RNA <- rbind(PARNA1_RNA_relative, PARNA2_RNA_relative, PARNA3_RNA_relative, PARNA4_RNA_relative)
#Substitute the name for BBBB with Control
PARNA_RNA$RE_Identity <- gsub("BBBB", "Control", PARNA_RNA$RE_Identity)
#make a new column called "condensed" with all "B" removed
PARNA_RNA$condensed <- gsub("B", "", PARNA_RNA$RE_Identity)
#group by RE_Identity and the summarise
PARNA_RNA_grouped <- dplyr::group_by(PARNA_RNA, RE_Identity) 
PARNA_RNA_grouped$FOLD_linear <- 2^PARNA_RNA_grouped$FOLD
PARNA_RNA_relative_sum <- dplyr::summarise(PARNA_RNA_grouped, Fold_TE=median(FOLD), mad_TE=mad(FOLD, constant=1), Fold_TE_mean=mean(FOLD), sd_TE=sd(FOLD), se_TE=se(FOLD), norm_mean_TE=mean(Normalized), sd_norm_TE=sd(Normalized), lin_mean_TE=mean(FOLD_linear), sd_lin_TE=sd(FOLD_linear), iqr_TE=IQR(FOLD), count_TE=length(FOLD)/4, q1_TE=unname(quantile(FOLD,0.25)), q3_TE=unname(quantile(FOLD, 0.75)))

#group by Condensed and the summarise
PARNA_RNA_grouped_2 <- dplyr::group_by(PARNA_RNA, condensed) 
PARNA_RNA_condensed_sum <- dplyr::summarise(PARNA_RNA_grouped_2, Fold_TE=median(FOLD), mad_TE=mad(FOLD, constant=1), Fold_TE_mean=mean(FOLD), sd_TE=sd(FOLD))

#find replicate results for PARNA normalized by plasmid
PARNA_plasmid <- rbind(PARNA1_plasmid_relative, PARNA2_plasmid_relative, PARNA3_plasmid_relative, PARNA4_plasmid_relative)
#Substitute the name for BBBB with Control
PARNA_plasmid $RE_Identity <- gsub("BBBB", "Control", PARNA_plasmid $RE_Identity)
#make a new column called "condensed" with all "B" removed
PARNA_plasmid $condensed <- gsub("B", "", PARNA_plasmid $RE_Identity)
#group by RE_Identity and the summarise
PARNA_plasmid_grouped <- dplyr::group_by(PARNA_plasmid, RE_Identity) 
PARNA_plasmid_relative_sum <- dplyr::summarise(PARNA_plasmid_grouped, Fold_RNA_PS=median(FOLD), mad_RNA_PS=mad(FOLD, constant=1), Fold_RNA_PS_mean=mean(FOLD), sd_RNA_PS=sd(FOLD))

#group by Condensed and the summarise
PARNA_plasmid_grouped_2 <- dplyr::group_by(PARNA_plasmid, condensed) 
PARNA_plasmid_condensed_sum <- dplyr::summarise(PARNA_plasmid_grouped_2, Fold_RNA_PS=median(FOLD), mad_RNA_PS=mad(FOLD, constant=1), Fold_RNA_PS_mean=mean(FOLD), sd_RNA_PS=sd(FOLD))

#bind replicate results for RNA
RNA_plasmid <- rbind(RNA1_plasmid_relative, RNA2_plasmid_relative, RNA3_plasmid_relative, RNA4_plasmid_relative)
#Substitute the name for BBBB with Control
RNA_plasmid$RE_Identity <- gsub("BBBB", "Control", RNA_plasmid$RE_Identity)
#make a new column called "condensed" with all "B" removed
RNA_plasmid$condensed <- gsub("B", "", RNA_plasmid$RE_Identity)
#group by RE_Identity and the summarise
RNA_plasmid_grouped <- dplyr::group_by(RNA_plasmid, RE_Identity) 
RNA_plasmid_grouped$FOLD_linear <- 2^RNA_plasmid_grouped$FOLD
RNA_plasmid_relative_sum <- dplyr::summarise(RNA_plasmid_grouped, Fold_RNA=median(FOLD), mad_RNA=mad(FOLD, constant=1), Fold_RNA_mean=mean(FOLD), sd_RNA=sd(FOLD), se_RNA=se(FOLD), norm_mean_RNA=mean(Normalized), sd_norm_RNA=sd(Normalized), lin_mean_RNA=mean(FOLD_linear), sd_lin_RNA=sd(FOLD_linear), iqr_RNA=IQR(FOLD), count_RNA=length(FOLD)/4, q1_RNA=unname(quantile(FOLD,0.25)), q3_RNA=unname(quantile(FOLD, 0.75)))

#group by Condensed and the summarise
RNA_plasmid_grouped_2 <- dplyr::group_by(RNA_plasmid, condensed) 
RNA_plasmid_condensed_sum <- dplyr::summarise(RNA_plasmid_grouped_2, Fold_RNA=median(FOLD), mad_RNA=mad(FOLD, constant=1), Fold_RNA_mean=mean(FOLD), sd_RNA=sd(FOLD))


#merge TE and RNA results
TE_RNA <- merge(PARNA_RNA_relative_sum, RNA_plasmid_relative_sum, by.x="RE_Identity", by.y="RE_Identity")
#remove Let-7 binding site variants to make data.frame with only 625 combinations
TE_RNA_subset <- subset(TE_RNA, !RE_Identity %in% c("6merx1", "6merx2", "6merx4", "7mera1x1", "7mera1x2", "7mera1x4" ,"7merm8x1", "7merm8x2", "7merm8x4", "7pcx2", "7pcx4", "7pcx1", "Dna2", "FIGNL2", "SMARCAD1", "HMGA2", "C14orf28"))
#rename control
TE_RNA_subset$RE_Identity <- gsub("Control", "BBBB", TE_RNA_subset$RE_Identity)
#subset RNA results or TE results
RNA_subset <- data.frame(RE_Identity = TE_RNA_subset$RE_Identity, Fold = TE_RNA_subset$Fold_RNA)
TE_subset <- data.frame(RE_Identity = TE_RNA_subset$RE_Identity, Fold = TE_RNA_subset$Fold_TE)
#write tables of summarized RNA and TE results
write.table(RNA_subset, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/RNA_summarized.txt")
write.table(TE_subset, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/TE_summarized.txt")
#add PARNA normalized by plasmid results to TE and RNA results
TE_RNA_PS <- merge(TE_RNA, PARNA_plasmid_relative_sum, by.x="RE_Identity", by.y="RE_Identity")
TE_RNA_relative_sum <- TE_RNA_PS
#make new column called group based on RE_Identity
TE_RNA_relative_sum$group <- TE_RNA_relative_sum$RE_Identity
#Remove "B" from group
TE_RNA_relative_sum$group <- gsub("B", "", TE_RNA_relative_sum$group)
#write table of TE and RNA summarized data
write.table(TE_RNA_relative_sum, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/TE_RNA.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
#load TE_RNA.txt file that has been annotated by hand
annotated <- read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/TE_RNA_annotated.txt")
#merge annotated data with TE_RNA_relative_sum
TE_RNA_annotated <- merge(TE_RNA_relative_sum, annotated, by="RE_Identity")
#load cbpallete and theme_science
cbPalette <- c("#D55E00", "gray", "black", "#56B4E9", "#009E73", "#E69F00")
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

#Rename values in Group
TE_RNA_annotated$group.y <- gsub("HuR", "ARE", TE_RNA_annotated$group.y)
TE_RNA_annotated$group.y <- gsub("Pumilio", "PRE", TE_RNA_annotated$group.y)
TE_RNA_annotated$group.y <- gsub("Smaug", "SRE", TE_RNA_annotated$group.y)
#define legend order
legend_order <- c("Control", "Combination", "Let-7", "PRE", "ARE", "SRE")
#plot of TE vs RNA fold change
ggplot(TE_RNA_annotated, aes(x = Fold_RNA, y= Fold_TE, colour=group.y)) + 
  scale_colour_manual(name = "CRE", values=cbPalette, breaks = legend_order) + 
  xlab("Relative Expression (log2)") + ylab("Relative TE (log2)") +
  geom_hline(yintercept = 0) +  geom_vline(xintercept = 0) +
  geom_point(alpha=0.5) + theme_science()
ggsave("TE_RNA.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)
#make new data.frame for TE summarized data
PARNA_RNA_relative_sum_rename <- PARNA_RNA_relative_sum
PARNA_RNA_relative_sum_rename$Group <- "TE"
PARNA_RNA_relative_sum_rename <- data.frame(RE_Identity = PARNA_RNA_relative_sum_rename$RE_Identity, Group = PARNA_RNA_relative_sum_rename$Group, Fold= PARNA_RNA_relative_sum_rename$Fold_TE)
#make new data.frame for RNA summarized data
RNA_plasmid_relative_sum_rename <- RNA_plasmid_relative_sum
RNA_plasmid_relative_sum_rename$Group <- "RNA"
RNA_plasmid_relative_sum_rename <- data.frame(RE_Identity = RNA_plasmid_relative_sum_rename$RE_Identity, Group = RNA_plasmid_relative_sum_rename$Group, Fold= RNA_plasmid_relative_sum_rename$Fold_RNA)
#bind TE and RNA summarized data and write table
TE_RNA_bound <- rbind(PARNA_RNA_relative_sum_rename, RNA_plasmid_relative_sum_rename)
write.table(TE_RNA_bound, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/TE_RNA_bound.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
#plot TE and RNA fold change 
ggplot(TE_RNA_bound, aes(x=reorder(RE_Identity, -Fold), y=Fold, colour=Group)) + 
  scale_colour_manual(values= c("grey", "black"), name="Level") +
  geom_point() + xlab("RE Identity") + 
  ylab("Fold Change (log2)") + theme_science() + 
  geom_hline(yintercept = 0, linetype="dashed") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("Boom.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)


#merge RNA and TE by condensed and then merge wiht PARNA normalized by plasmid
TE_RNA_condensed <- merge(PARNA_RNA_condensed_sum, RNA_plasmid_condensed_sum, by.x="condensed", by.y="condensed")
TE_RNA_PS_condensed <- merge(TE_RNA_condensed, PARNA_plasmid_condensed_sum, by.x="condensed", by.y="condensed")

write.table(TE_RNA_PS_condensed, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Results_Condensed.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
#write table of RNA, TE and PARNA normalized plasmid unsummarized data
write.table(RNA_plasmid, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/RNA_plasmid.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(PARNA_plasmid, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/PARNA_plasmid.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(PARNA_RNA, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/PARNA_RNA.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)


#Open "Grouped_barcode_Identities"
Grouped_barcode_identities = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Grouped_Barcode_identities.txt")
#Filter out RE containing only Pumilio and Blank, make barplot of Relative Expression and write txt file
pumilio <- dplyr::filter(TE_RNA_relative_sum, RE_Identity %in% Grouped_barcode_identities$Pumilio)
write.table(pumilio, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Pumilio_Results.txt", sep="\t")
#Filter out RE containing only Let-7 and Blank, make barplot of Relative Expression and write txt file
let_7 <- dplyr::filter(TE_RNA_relative_sum, RE_Identity %in% Grouped_barcode_identities$Let.7)
write.table(let_7, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Let-7_Results.txt",  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE) 
#Filter out RE containing only Natural targets and Blank, make barplot of Relative Expression and write txt file
Natural_targets <- dplyr::filter(TE_RNA_relative_sum, RE_Identity %in% Grouped_barcode_identities$Natural)
write.table(Natural_targets, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Natural_Targets_Results.txt",  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
#Filter out RE containing only Smaug and Blank, make barplot of Relative Expression and write txt file
smaug <- dplyr::filter(TE_RNA_relative_sum, RE_Identity %in% Grouped_barcode_identities$Smg)
write.table(smaug, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Smaug_Results.txt",  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
#Filter out RE containing only HuR and Blank, make barplot of Relative Expression and write txt file
HuR <- dplyr::filter(TE_RNA_relative_sum, RE_Identity %in% Grouped_barcode_identities$HuR)
write.table(HuR, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/HuR_Results.txt",  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
#Filter out RE containing only Altered Seed sequence and Blank, make barplot of Relative Expression and write txt file
Seed <- dplyr::filter(TE_RNA_relative_sum, RE_Identity %in% Grouped_barcode_identities$Seed)
write.table(Seed, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Seed_Results.txt",  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
#Filter out RE containing Pum, Let-7 and Blank sites
Pum_Let_7 <- dplyr::filter(TE_RNA_relative_sum, RE_Identity %in% Grouped_barcode_identities$Pum_Let_7)
write.table(Pum_Let_7, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Pum_Let_7_Results.txt",  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

#Filter out RE containing HuR, Let-7 and Blank sites
HuR_Let_7 <- dplyr::filter(TE_RNA_relative_sum, RE_Identity %in% Grouped_barcode_identities$HuR_Let_7)
write.table(HuR_Let_7, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/HuR_Let_7_Results.txt",  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

#Filter out RE containing Smg, Let-7 and blank sites
Smg_Let_7 <- dplyr::filter(TE_RNA_relative_sum, RE_Identity %in% Grouped_barcode_identities$Smg_Let_7)
write.table(Smg_Let_7, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Smg_Let_7_Results.txt",  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)


#Filter out RE containing HuR, Pum and Blank sites
HuR_Pum <- dplyr::filter(TE_RNA_relative_sum, RE_Identity %in% Grouped_barcode_identities$HuR_Pum)
write.table(HuR_Pum, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/HuR_Pum_Results.txt",  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

#Filter out Smaug, Pum and Blank sites
Smg_Pum <- subset(TE_RNA_relative_sum, grepl("Control|[SBp][SBp][SBp][SBp]", TE_RNA_relative_sum$RE_Identity))
write.table(Smg_Pum, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Smg_Pum_Results.txt",  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

#Filter out Smaug, Pum and Blank sites
Let_PC <- dplyr::filter(TE_RNA_relative_sum, RE_Identity %in% Grouped_barcode_identities$Let_PC)
write.table(Smg_Pum, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Let_PC_Results.txt",  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)

#calculate coefficient of variance for RNA expression, must use mean fold change (linear not log)
RNA_plasmid_relative_sum$CV <- abs(RNA_plasmid_relative_sum$sd_lin_RNA/RNA_plasmid_relative_sum$lin_mean_RNA)
#determine mean and median CV and plot density and histogram of CV values
mean(RNA_plasmid_relative_sum$CV)
median(RNA_plasmid_relative_sum$CV)
ggplot(RNA_plasmid_relative_sum, aes(CV)) + geom_density() + theme_science()
ggplot(RNA_plasmid_relative_sum, aes(CV)) + geom_histogram() + ylab("Count")+ theme_science() + scale_y_discrete(limits=(0:150), breaks=c(0,50,100, 150))
ggsave("RNA_CV.tiff", width=5, height=3, units="in")

#determine mean and median SD for RNA and plot density and histogram
mean(RNA_plasmid_relative_sum$sd_RNA)
median(RNA_plasmid_relative_sum$sd_RNA)
ggplot(RNA_plasmid_relative_sum, aes(sd_RNA)) + geom_density()
ggplot(RNA_plasmid_relative_sum, aes(sd_RNA)) + geom_histogram()
#determine mean and median SE for RNA and plot density and histogram
mean(RNA_plasmid_relative_sum$se_RNA)
median(RNA_plasmid_relative_sum$se_RNA)
ggplot(RNA_plasmid_relative_sum, aes(se_RNA)) + geom_density()
ggplot(RNA_plasmid_relative_sum, aes(se_RNA)) + geom_histogram()
#determine mean and median IQR for RNA and plot density and histogram
mean(RNA_plasmid_relative_sum$iqr_RNA)
median(RNA_plasmid_relative_sum$iqr_RNA)
ggplot(RNA_plasmid_relative_sum, aes(iqr_RNA)) + geom_density()
ggplot(RNA_plasmid_relative_sum, aes(iqr_RNA)) + geom_histogram()
#determine mean and median for fold change of RNA and plot density and histogram
mean(RNA_plasmid_relative_sum$Fold_RNA)
median(RNA_plasmid_relative_sum$Fold_RNA)
ggplot(RNA_plasmid_relative_sum, aes(Fold_RNA)) + geom_density()
ggplot(RNA_plasmid_relative_sum, aes(Fold_RNA)) + geom_histogram()

#subset data to remove barcodes counted >10 (control and let-7-pcx4 have 50 barcodes each)
RNA_plasmid_relative_sum_counts <- subset(RNA_plasmid_relative_sum, !count_RNA >10)
#determine mean and median for number of barcodes counted and plot density and histogram
mean(RNA_plasmid_relative_sum_counts$count_RNA)
median(RNA_plasmid_relative_sum_counts$count_RNA)
ggplot(RNA_plasmid_relative_sum_counts, aes(count_RNA)) + geom_density()
ggplot(RNA_plasmid_relative_sum_counts, aes(count_RNA)) + geom_histogram()

#calculate coefficient of variance for RNA expression, must use mean fold change (linear not log)
PARNA_RNA_relative_sum$CV <- abs(PARNA_RNA_relative_sum$sd_lin_TE/PARNA_RNA_relative_sum$lin_mean_TE)
#determine mean and median CV and plot density and histogram of CV values
mean(PARNA_RNA_relative_sum$CV)
median(PARNA_RNA_relative_sum$CV)
ggplot(PARNA_RNA_relative_sum, aes(CV)) + geom_density()
ggplot(PARNA_RNA_relative_sum, aes(CV)) + geom_histogram() + ylab("Count")+ theme_science() + scale_y_discrete(limits=(0:150), breaks=c(0,20,40,60,80))
ggsave("TE_CV.tiff", width=5, height=3, units="in")

#determine mean and median SD, SE and Fold change for TE
mean(PARNA_RNA_relative_sum$sd_TE)
median(PARNA_RNA_relative_sum$sd_TE)
mean(PARNA_RNA_relative_sum$se_TE)
median(PARNA_RNA_relative_sum$se_TE)
mean(PARNA_RNA_relative_sum$Fold_TE)
median(PARNA_RNA_relative_sum$Fold_TE)
mean(PARNA_RNA_relative_sum$Fold_TE)
median(PARNA_RNA_relative_sum$Fold_TE)

#subset data to remove barcodes counted >10 (control and let-7-pcx4 have 50 barcodes each)
PARNA_RNA_relative_sum_counts <- subset(PARNA_RNA_relative_sum, !count_TE >10)
#determine mean and median for number of barcodes counted and plot histogram
mean(PARNA_RNA_relative_sum$count_TE)
median(PARNA_RNA_relative_sum$count_TE)
ggplot(PARNA_RNA_relative_sum_counts, aes(count_TE)) + geom_histogram()

#merge RNA_plasmid grouped (unsummarized) with RNA plasmid summarized to get FOLD change of each barcode or replicated along with iqr and quartile values for each RE
RNA_plasmid_grouped_outlier <- merge(RNA_plasmid_grouped, RNA_plasmid_relative_sum, by="RE_Identity")
#drop outliers
RNA_plasmid_grouped_outlier <-RNA_plasmid_grouped_outlier[!(RNA_plasmid_grouped_outlier$FOLD < RNA_plasmid_grouped_outlier$q1_RNA-1.5*RNA_plasmid_grouped_outlier$iqr_RNA | RNA_plasmid_grouped_outlier$FOLD > RNA_plasmid_grouped_outlier$q3_RNA+1.5*RNA_plasmid_grouped_outlier$iqr_RNA),]
#group by RE_Identity and summarize values
RNA_plasmid_grouped_outlier <- dplyr::group_by(RNA_plasmid_grouped_outlier, RE_Identity)
RNA_plasmid_outlier_sum <- dplyr::summarise(RNA_plasmid_grouped_outlier, lin_mean_RNA=mean(FOLD_linear), sd_lin_RNA=sd(FOLD_linear), count_RNA=length(FOLD)/4, Fold = median(FOLD), SD=sd(FOLD), mean_Fold= mean(FOLD), SE=se(FOLD))

#calculate coefficient of variance for RNA expression, must use mean fold change (linear not log)
RNA_plasmid_outlier_sum$CV <- abs(RNA_plasmid_outlier_sum$sd_lin_RNA/RNA_plasmid_outlier_sum$lin_mean_RNA)
#find mean and median CV and make plots
mean(RNA_plasmid_outlier_sum$CV)
median(RNA_plasmid_outlier_sum$CV)
ggplot(RNA_plasmid_outlier_sum, aes(CV)) + geom_density()
ggplot(RNA_plasmid_outlier_sum, aes(CV)) + geom_histogram() + ylab("Count")+ theme_science() + scale_y_discrete(limits=(0:150), breaks=c(0,20,40,60,80))
ggsave("RNA_CV_no_outliers.tiff", width=5, height=3, units="in")

#subset data to remove barcodes counted >10 (control and let-7-pcx4 have 50 barcodes each)
RNA_plasmid_outlier_sum_counts <- subset(RNA_plasmid_outlier_sum, !count_RNA >10)
#find mean and median barcodes counted and make plots
mean(RNA_plasmid_outlier_sum_counts$count_RNA)
median(RNA_plasmid_outlier_sum_counts$count_RNA)
ggplot(RNA_plasmid_outlier_sum_counts, aes(count_RNA)) + geom_density()
ggplot(RNA_plasmid_outlier_sum_counts, aes(count_RNA)) + geom_histogram()

#compare data with and without outliers
RNA_sum_total <- merge(RNA_plasmid_relative_sum, RNA_plasmid_outlier_sum, by="RE_Identity")
ggplot(RNA_sum_total, aes(Fold_RNA, Fold)) + geom_point() + geom_smooth(method="lm")
cor.test(RNA_sum_total$Fold_RNA, RNA_sum_total$Fold)

#determine mean and median SD, SE and Fold change for RNA
mean(RNA_plasmid_outlier_sum_counts$mean_Fold)
mean(RNA_plasmid_outlier_sum_counts$SD)
mean(RNA_plasmid_outlier_sum_counts$Fold)
mean(RNA_plasmid_outlier_sum_counts$SE)
mean(RNA_plasmid_outlier_sum_counts$count_RNA)
median(RNA_plasmid_outlier_sum_counts$mean_Fold)
median(RNA_plasmid_outlier_sum_counts$SD)
median(RNA_plasmid_outlier_sum_counts$Fold)
median(RNA_plasmid_outlier_sum_counts$SE)
median(RNA_plasmid_outlier_sum_counts$count_RNA)

#merge RNA_plasmid grouped (unsummarized) with RNA plasmid summarized to get FOLD change of each barcode or replicated along with iqr and quartile values for each RE
PARNA_RNA_grouped_outlier <- merge(PARNA_RNA_grouped, PARNA_RNA_relative_sum, by="RE_Identity")
PARNA_RNA_grouped_outlier <-PARNA_RNA_grouped_outlier[!(PARNA_RNA_grouped_outlier$FOLD < PARNA_RNA_grouped_outlier$q1_TE-1.5*PARNA_RNA_grouped_outlier$iqr_TE | PARNA_RNA_grouped_outlier$FOLD > PARNA_RNA_grouped_outlier$q3_TE+1.5*PARNA_RNA_grouped_outlier$iqr_TE),]
PARNA_RNA_grouped_outlier <- dplyr::group_by(PARNA_RNA_grouped_outlier, RE_Identity)
PARNA_RNA_outlier_sum <- dplyr::summarise(PARNA_RNA_grouped_outlier, lin_mean_RNA=mean(FOLD_linear), sd_lin_RNA=sd(FOLD_linear), count_RNA=length(FOLD)/4, Fold = median(FOLD), SD=sd(FOLD), mean_Fold= mean(FOLD), SE=se(FOLD)) 

#calculate coefficient of variance for TE, must use mean fold change (linear not log)
PARNA_RNA_outlier_sum$CV <- abs(PARNA_RNA_outlier_sum$sd_lin_RNA/PARNA_RNA_outlier_sum$lin_mean_RNA)
mean(PARNA_RNA_outlier_sum$CV)
median(PARNA_RNA_outlier_sum$CV)
ggplot(PARNA_RNA_outlier_sum, aes(CV)) + geom_density()
ggplot(PARNA_RNA_outlier_sum, aes(CV)) + geom_histogram() + ylab("Count")+ theme_science() + scale_y_discrete(limits=(0:150), breaks=c(0,10,20,30,40,50))
ggsave("TE_CV_no_outliers.tiff", width=5, height=3, units="in")

#subset data to remove barcodes counted >10 (control and let-7-pcx4 have 50 barcodes each)
PARNA_RNA_outlier_sum_counts <- subset(PARNA_RNA_outlier_sum, !count_RNA >10)
#find mean and median barcodes counted and make plots
mean(PARNA_RNA_outlier_sum_counts$count_RNA)
median(PARNA_RNA_outlier_sum_counts$count_RNA)
ggplot(PARNA_RNA_outlier_sum_counts, aes(count_RNA)) + geom_density()
ggplot(PARNA_RNA_outlier_sum_counts, aes(count_RNA)) + geom_histogram()

#determine mean and median SD, SE and Fold change for TE
mean(PARNA_RNA_outlier_sum_counts$mean_Fold)
mean(PARNA_RNA_outlier_sum_counts$SD)
mean(PARNA_RNA_outlier_sum_counts$Fold)
mean(PARNA_RNA_outlier_sum_counts$SE)
mean(PARNA_RNA_outlier_sum_counts$count_RNA)
median(PARNA_RNA_outlier_sum_counts$mean_Fold)
median(PARNA_RNA_outlier_sum_counts$SD)
median(PARNA_RNA_outlier_sum_counts$Fold)
median(PARNA_RNA_outlier_sum_counts$SE)
median(PARNA_RNA_outlier_sum_counts$count_RNA)

#compare data with and without outliers
PARNA_RNA_sum_total <- merge(PARNA_RNA_relative_sum, PARNA_RNA_outlier_sum, by="RE_Identity")
ggplot(PARNA_RNA_sum_total, aes(Fold_TE, Fold)) + geom_point() + geom_smooth(method="lm")
cor.test(PARNA_RNA_sum_total$Fold_TE, PARNA_RNA_sum_total$Fold)

