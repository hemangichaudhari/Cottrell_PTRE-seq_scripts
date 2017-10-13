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



#open sequencing results
Plasmid_1_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L001_05_counts")
Plasmid_2_L1 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L001_05_counts")

Plasmid_1_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L002_05_counts")
Plasmid_2_L2 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L002_05_counts")

Plasmid_1_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L003_05_counts")
Plasmid_2_L3 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L003_05_counts")

Plasmid_1_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S18_L004_05_counts")
Plasmid_2_L4 = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/S17_L004_05_counts")

#open sequencing results from cell-lines
HEK = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Cells/Results_01_counts")
HDF = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Cells/Results_02_counts")
N2A = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Cells/Results_03_counts")

#bind plasmid sequencing results
Plasmid_1 = rbind(Plasmid_1_L1, Plasmid_1_L2, Plasmid_1_L3, Plasmid_1_L4)
Plasmid_1 <- dplyr::group_by(Plasmid_1, V2)
Plasmid_1 <- dplyr::summarise(Plasmid_1, Counts=sum(V1))

Plasmid_2 = rbind(Plasmid_2_L1, Plasmid_2_L2, Plasmid_2_L3, Plasmid_2_L4)
Plasmid_2 <- dplyr::group_by(Plasmid_2, V2)
Plasmid_2 <- dplyr::summarise(Plasmid_2, Counts=sum(V1))

#change column names
colnames(Plasmid_1) = c("Barcode", "Counts_Plasmid")
colnames(Plasmid_2) = c("Barcode", "Counts_Plasmid")
colnames(HEK) = c("Counts_cDNA", "Barcode")
colnames(HDF) = c("Counts_cDNA", "Barcode")
colnames(N2A) = c("Counts_cDNA", "Barcode")

#open "Barcode_Identities"
barcodes = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Cells/Barcode_identities.txt")

#join "Barcode_Identities" with Sequencing Results
Plasmid_1_id = merge(Plasmid_1, barcodes, by="Barcode")
Plasmid_2_id = merge(Plasmid_2, barcodes, by="Barcode")

HEK_id = merge(HEK,barcodes, by ="Barcode")
HDF_id = merge(HDF,barcodes, by ="Barcode")
N2A_id = merge(N2A,barcodes, by ="Barcode")


#Omit barcodes with no reads (NA)
Plasmid_1_id = na.omit(Plasmid_1_id)
Plasmid_2_id = na.omit(Plasmid_2_id)

#Omit barcodes with no reads (NA)
HEK_id = na.omit(HEK_id)

HEK_all <- HEK_id

#rename some of the regulatory elements
HEK_all$RE_Identity <- gsub("h", "A", HEK_all$RE_Identity)
HEK_all$RE_Identity <- gsub("7", "L", HEK_all$RE_Identity)
HEK_all$RE_Identity <- gsub("Lmer", "7mer", HEK_all$RE_Identity)

#Omit barcodes with no reads (NA)
HDF_id = na.omit(HDF_id)

HDF_all <- HDF_id

#rename some of the regulatory elements
HDF_all$RE_Identity <- gsub("h", "A", HDF_all$RE_Identity)
HDF_all$RE_Identity <- gsub("7", "L", HDF_all$RE_Identity)
HDF_all$RE_Identity <- gsub("Lmer", "7mer", HDF_all$RE_Identity)

#Omit barcodes with no reads (NA)
N2A_id = na.omit(N2A_id)

N2A_all <- N2A_id
#rename some of the regulatory elements
N2A_all$RE_Identity <- gsub("h", "A", N2A_all$RE_Identity)
N2A_all$RE_Identity <- gsub("7", "L", N2A_all$RE_Identity)
N2A_all$RE_Identity <- gsub("Lmer", "7mer", N2A_all$RE_Identity)

#assign cellt ype to each data set
HEK_all$Cell_type <- "HEK"
HDF_all$Cell_type <- "HDF"
N2A_all$Cell_type <- "N2A"

#bind sequencing results for each cell type and write table
cells_all <- rbind(HEK_all, HDF_all, N2A_all)

write.table(cells_all , "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/cells_reads.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)


#Omit barcodes with less than 10 reads not functional at this time

Plasmid_1_id = subset(Plasmid_1_id, Counts_Plasmid > 10)
Plasmid_2_id = subset(Plasmid_2_id, Counts_Plasmid > 10)

#merge plasmid sequencing replicates
Plasmid <- merge(Plasmid_1_id, Plasmid_2_id, by.x = "Barcode", by.y = "Barcode")

#determine average couts for each plasmid
Plasmid$Counts_avg <- (Plasmid$Counts_Plasmid.x + Plasmid$Counts_Plasmid.y)/2

#Group barcodes by RE_Identity
Plasmid  <- dplyr::group_by(Plasmid , RE_Identity.x) 

HEK_id <- dplyr::group_by(HEK_id, RE_Identity)
HDF_id <- dplyr::group_by(HDF_id, RE_Identity)
N2A_id <- dplyr::group_by(N2A_id, RE_Identity)


#join cDNA and plasmid by barcode then omit "na"

HEK_plasmid = merge(HEK_id, Plasmid, by ="Barcode")
HEK_plasmid_na = na.omit(HEK_plasmid)
HDF_plasmid = merge(HDF_id,Plasmid, by ="Barcode")
HDF_plasmid_na = na.omit(HDF_plasmid)
N2A_plasmid = merge(N2A_id,Plasmid, by ="Barcode")
N2A_plasmid_na = na.omit(N2A_plasmid)


#normalize cDNA_counts to plasmid_counts
HEK_plasmid_na$normalized = HEK_plasmid_na$Counts_cDNA/HEK_plasmid_na$Counts_avg
HDF_plasmid_na$normalized = HDF_plasmid_na$Counts_cDNA/HDF_plasmid_na$Counts_avg
N2A_plasmid_na$normalized = N2A_plasmid_na$Counts_cDNA/N2A_plasmid_na$Counts_avg


#make new data.frame with just barcode identities, and normalized values
HEK_plasmid_normalized = data.frame(HEK_plasmid_na$RE_Identity.x, HEK_plasmid_na$normalized)
colnames(HEK_plasmid_normalized) = c("RE_Identity","Normalized")
HDF_plasmid_normalized = data.frame(HDF_plasmid_na$RE_Identity.x, HDF_plasmid_na$normalized)
colnames(HDF_plasmid_normalized) = c("RE_Identity","Normalized")
N2A_plasmid_normalized = data.frame(N2A_plasmid_na$RE_Identity.x, N2A_plasmid_na$normalized)
colnames(N2A_plasmid_normalized) = c("RE_Identity","Normalized")


#Define 'x' as the value for the control "BBBB"
x1 = subset(HEK_plasmid_normalized, RE_Identity=="BBBB")
x2 = subset(HDF_plasmid_normalized, RE_Identity=="BBBB")
x3 = subset(N2A_plasmid_normalized, RE_Identity=="BBBB")

#Add a new collumn with relative data for all RE and name it "Fold"
HEK_plasmid_relative <- dplyr::mutate(HEK_plasmid_normalized, log2(Normalized/median(x1$Normalized)))
colnames(HEK_plasmid_relative)[colnames(HEK_plasmid_relative)=="log2(Normalized/median(x1$Normalized))"] <- "FOLD"
HEK_plasmid_control <- HEK_plasmid_relative
#rename control reporter
HEK_plasmid_control$RE_Identity <- gsub("BBBB", "Control", HEK_plasmid_control$RE_Identity)

#repeat as above for HDF
HDF_plasmid_relative <- dplyr::mutate(HDF_plasmid_normalized, log2(Normalized/median(x2$Normalized)))
colnames(HDF_plasmid_relative)[colnames(HDF_plasmid_relative)=="log2(Normalized/median(x2$Normalized))"] <- "FOLD"
HDF_plasmid_control <- HDF_plasmid_relative
HDF_plasmid_control$RE_Identity <- gsub("BBBB", "Control", HDF_plasmid_control$RE_Identity)

#repeat as above for N2A
N2A_plasmid_relative <- dplyr::mutate(N2A_plasmid_normalized, log2(Normalized/median(x3$Normalized)))
colnames(N2A_plasmid_relative)[colnames(N2A_plasmid_relative)=="log2(Normalized/median(x3$Normalized))"] <- "FOLD"
N2A_plasmid_control <- N2A_plasmid_relative
N2A_plasmid_control$RE_Identity <- gsub("BBBB", "Control", N2A_plasmid_control$RE_Identity)

#group by RE_Identity and the summarise
HEK_plasmid_relative_grouped <- dplyr::group_by(HEK_plasmid_relative, RE_Identity) 
HEK_plasmid_relative_sum <- dplyr::summarise(HEK_plasmid_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1))

HDF_plasmid_relative_grouped <- dplyr::group_by(HDF_plasmid_relative, RE_Identity) 
HDF_plasmid_relative_sum <- dplyr::summarise(HDF_plasmid_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1))

N2A_plasmid_relative_grouped <- dplyr::group_by(N2A_plasmid_relative, RE_Identity) 
N2A_plasmid_relative_sum <- dplyr::summarise(N2A_plasmid_relative_grouped, Fold=median(FOLD), max=max(FOLD), min=min(FOLD), mad=mad(FOLD, constant=1))

#rename B as *

HEK_plasmid_control$condensed <- gsub("B", "", HEK_plasmid_control$RE_Identity)
HDF_plasmid_control$condensed <- gsub("B", "", HDF_plasmid_control$RE_Identity)
N2A_plasmid_control$condensed <- gsub("B", "", N2A_plasmid_control$RE_Identity)

write.table(HEK_plasmid_control, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_Cells/HEK_plasmid.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(HDF_plasmid_control, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_Cells/HDF_plasmid.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(N2A_plasmid_control, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_Cells/N2A_plasmid.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)


