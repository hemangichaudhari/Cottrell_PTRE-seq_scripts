library(dplyr)
library(ggplot2)
library(extrafont)
library(stringr)
library(gridExtra)
library(grid)
font_import()

loadfonts(device="win")
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

#load RNA results
RNA_results= read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/RNA_plasmid.txt")
#load grouped_barcode identities
Grouped_barcode_identities = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Grouped_Barcode_identities.txt")
#subset RNA results for only those targeted by Let-7 and/or pumilio
RNA_results_pum_let <- dplyr::filter(RNA_results, RE_Identity %in% Grouped_barcode_identities$Pum_Let_7)
#subset Pum_let to compare position of 1 of either element
position_pum_let <- subset(RNA_results_pum_let, RE_Identity %in% c("7BBB", "7BBp", "7BpB", "7pBB", "pBBB", "pBB7", "pB7B", "p7BB"))
#rename 7-> L and B ->*
position_pum_let$RE_Identity <- gsub("7", "L", position_pum_let$RE_Identity)
position_pum_let$RE_Identity <- gsub("B", "*", position_pum_let$RE_Identity)
#plot position_pum_let
ggplot(position_pum_let, aes(RE_Identity, FOLD)) + 
  geom_boxplot(outlier.colour="white", notch = TRUE, fill="grey") + scale_y_continuous(limits=c(-3.8,0.5)) +
  ylab("Relative Expression (log2)") + xlab("RE Identity") + 
  geom_hline(yintercept = 0, linetype="dashed") + theme_science()
ggsave("Pum_Let_position_RNA.tiff", height=3, width=5, units="in")

#Make subsets of the the data for all arrangements of N Let-7 sites and N PREs
Pum1_Let1 <- subset(RNA_results_pum_let, condensed %in% c('7p', 'p7', 'p', '7'))
Pum1_Let1$condensed <- gsub("7p|p7", "Comb.", Pum1_Let1$condensed)
Pum2_Let1 <- subset(RNA_results_pum_let, condensed %in% c('7pp', 'pp7', 'p7p', 'pp', '7'))
Pum2_Let1$condensed <- gsub("7pp|pp7|p7p", "Comb.", Pum2_Let1$condensed)
Pum3_Let1 <- subset(RNA_results_pum_let, condensed %in% c('7ppp', 'ppp7', 'p7pp', 'pp7p', 'ppp', '7'))
Pum3_Let1$condensed <- gsub("7ppp|ppp7|p7pp|pp7p", "Comb.", Pum3_Let1$condensed)
Pum1_Let2 <- subset(RNA_results_pum_let, condensed %in% c('77p', 'p77', '7p7', '77', 'p'))
Pum1_Let2$condensed <- gsub("77p|p77|7p7", "Comb.", Pum1_Let2$condensed)
Pum2_Let2 <- subset(RNA_results_pum_let, condensed %in% c('77pp', '7p7p', '7pp7', 'p77p', 'p7p7', 'pp77', '77', 'pp'))
Pum2_Let2$condensed <- gsub('77pp|7p7p|7pp7|p77p|p7p7|pp77', "Comb.", Pum2_Let2$condensed)
Pum1_Let3 <- subset(RNA_results_pum_let, condensed %in% c('777p', 'p777', '7p77', '77p7', '777', 'p'))
Pum1_Let3$condensed <- gsub("777p|p777|77p7|7p77", "Comb.", Pum1_Let3$condensed)

#assign a group identity for each subseted data.frame
Pum1_Let1$Group <- "Pum1_Let1"
Pum2_Let1$Group <- "Pum2_Let1"
Pum3_Let1$Group <- "Pum3_Let1"
Pum1_Let2$Group <- "Pum1_Let2"
Pum2_Let2$Group <- "Pum2_Let2"
Pum1_Let3$Group <- "Pum1_Let3"

#bind subsets from above
pum_let <- rbind(Pum1_Let1, Pum2_Let1, Pum2_Let2, Pum3_Let1, Pum1_Let2, Pum2_Let2, Pum1_Let3)
pum_let_sum <- data.frame(Group = pum_let$Group, Fold = pum_let$FOLD, RE_Identity = pum_let$RE_Identity, condensed = pum_let$condensed)
#rename 7 -> L
pum_let_sum$condensed <- gsub("7", "L", pum_let_sum$condensed)
#plot RNA expression of each group from above 
group1 <- ggplot(subset(pum_let_sum, Group == "Pum1_Let1"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("p", "L", "Comb.")) + scale_y_continuous(limits=c(-4.7,0.5)) +
  geom_boxplot(outlier.colour="white", fill="grey", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 
group2 <- ggplot(subset(pum_let_sum, Group == "Pum2_Let1"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("pp", "L", "Comb.")) + scale_y_continuous(limits=c(-4.7,0.5)) +
  geom_boxplot(outlier.colour="white", fill="grey", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 
group3 <- ggplot(subset(pum_let_sum, Group == "Pum3_Let1"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("ppp", "L", "Comb.")) + scale_y_continuous(limits=c(-4.7,0.5)) +
  geom_boxplot(outlier.colour="white", fill="grey", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 
group4 <- ggplot(subset(pum_let_sum, Group == "Pum1_Let2"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("p", "LL", "Comb.")) + scale_y_continuous(limits=c(-4.7,0.5)) +
  geom_boxplot(outlier.colour="white", fill="grey", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 
group5 <- ggplot(subset(pum_let_sum, Group == "Pum2_Let2"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("pp", "LL", "Comb.")) + scale_y_continuous(limits=c(-4.7,0.5)) +
  geom_boxplot(outlier.colour="white", fill="grey", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 
group6 <- ggplot(subset(pum_let_sum, Group == "Pum1_Let3"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("p", "LLL", "Comb.")) + scale_y_continuous(limits=c(-4.7,0.5)) +
  geom_boxplot(outlier.colour="white", fill="grey", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 

#arrange above plots
g_rna <-grid.arrange(arrangeGrob(group1, group2, group3, group4, group5, group6, ncol=3, left=textGrob("Relative Expression", rot=90,gp =gpar(fontsize=12, fontfamily="Arial Black"), vjust=0.5)))
ggsave("Pum_Let_RNA.tiff", g_rna, height=4, width=5, dpi=300, units="in")


#Make subsets of the the data for all arrangements of N Let-7 sites and N PREs, summarise Fold
Pum1_Let1 <- subset(RNA_results_pum_let, condensed %in% c('7p', 'p7', 'p', '7'))
Pum1_Let1$condensed <- gsub("7p|p7", "Comb.", Pum1_Let1$condensed)
Pum1_Let1 <- dplyr::group_by(Pum1_Let1, condensed)
Pum1_Let1 <- dplyr::summarise(Pum1_Let1, Fold=median(FOLD), mad=mad(FOLD))
Pum2_Let1 <- subset(RNA_results_pum_let, condensed %in% c('7pp', 'pp7', 'p7p', 'pp', '7'))
Pum2_Let1$condensed <- gsub("7pp|pp7|p7p", "Comb.", Pum2_Let1$condensed)
Pum2_Let1 <- dplyr::group_by(Pum2_Let1, condensed)
Pum2_Let1 <- dplyr::summarise(Pum2_Let1, Fold=median(FOLD), mad=mad(FOLD))
Pum3_Let1 <- subset(RNA_results_pum_let, condensed %in% c('7ppp', 'ppp7', 'p7pp', 'pp7p', 'ppp', '7'))
Pum3_Let1$condensed <- gsub("7ppp|ppp7|p7pp|pp7p", "Comb.", Pum3_Let1$condensed)
Pum3_Let1 <- dplyr::group_by(Pum3_Let1, condensed)
Pum3_Let1 <- dplyr::summarise(Pum3_Let1, Fold=median(FOLD), mad=mad(FOLD))
Pum1_Let2 <- subset(RNA_results_pum_let, condensed %in% c('77p', 'p77', '7p7', '77', 'p'))
Pum1_Let2$condensed <- gsub("77p|p77|7p7", "Comb.", Pum1_Let2$condensed)
Pum1_Let2 <- dplyr::group_by(Pum1_Let2, condensed)
Pum1_Let2 <- dplyr::summarise(Pum1_Let2, Fold=median(FOLD), mad=mad(FOLD))
Pum2_Let2 <- subset(RNA_results_pum_let, condensed %in% c('77pp', '7p7p', '7pp7', 'p77p', 'p7p7', 'pp77', '77', 'pp'))
Pum2_Let2$condensed <- gsub('77pp|7p7p|7pp7|p77p|p7p7|pp77', "Comb.", Pum2_Let2$condensed)
Pum2_Let2 <- dplyr::group_by(Pum2_Let2, condensed)
Pum2_Let2 <- dplyr::summarise(Pum2_Let2, Fold=median(FOLD), mad=mad(FOLD))
Pum1_Let3 <- subset(RNA_results_pum_let, condensed %in% c('777p', 'p777', '7p77', '77p7', '777', 'p'))
Pum1_Let3$condensed <- gsub("777p|p777|77p7|7p77", "Comb.", Pum1_Let3$condensed)
Pum1_Let3 <- dplyr::group_by(Pum1_Let3, condensed)
Pum1_Let3 <- dplyr::summarise(Pum1_Let3, Fold=median(FOLD), mad=mad(FOLD))

#assign a group identity for each subseted data.frame
Pum1_Let1$Group <- "Pum1_Let1"
Pum2_Let1$Group <- "Pum2_Let1"
Pum3_Let1$Group <- "Pum3_Let1"
Pum1_Let2$Group <- "Pum1_Let2"
Pum2_Let2$Group <- "Pum2_Let2"
Pum1_Let3$Group <- "Pum1_Let3"

#bind subsets from above
pum_let <- rbind(Pum1_Let1, Pum2_Let1, Pum2_Let2, Pum3_Let1, Pum1_Let2, Pum2_Let2, Pum1_Let3)
pum_let_sum <- data.frame(Group = pum_let$Group, Fold = pum_let$Fold, condensed = pum_let$condensed)
#bind subsets from above
pum_let_sum$condensed <- gsub("7", "L", pum_let_sum$condensed)

#determine predicted expressin (sum of FOLD change of Let-7 and Pumilio)
Test <- "Predicted"
Pum1_Let1_ex <- Pum1_Let1$Fold[Pum1_Let1$condensed == "7"]+Pum1_Let1$Fold[Pum1_Let1$condensed == "p"]
Pum2_Let1_ex <- Pum2_Let1$Fold[Pum2_Let1$condensed == "7"]+Pum2_Let1$Fold[Pum2_Let1$condensed == "pp"]
Pum3_Let1_ex <- Pum3_Let1$Fold[Pum3_Let1$condensed == "7"]+Pum3_Let1$Fold[Pum3_Let1$condensed == "ppp"]
Pum1_Let2_ex <- Pum1_Let2$Fold[Pum1_Let2$condensed == "77"]+Pum1_Let2$Fold[Pum1_Let2$condensed == "p"]
Pum2_Let2_ex <- Pum2_Let2$Fold[Pum2_Let2$condensed == "77"]+Pum2_Let2$Fold[Pum2_Let2$condensed == "pp"]
Pum1_Let3_ex <- Pum1_Let3$Fold[Pum1_Let3$condensed == "777"]+Pum1_Let3$Fold[Pum1_Let3$condensed == "p"]
#make a list of above values
Fold <- c(Pum1_Let1_ex, Pum2_Let1_ex, Pum3_Let1_ex, Pum1_Let2_ex, Pum2_Let2_ex, Pum1_Let3_ex)
#make a list of group identities
Group <- c("Pum1_Let1", "Pum2_Let1", "Pum3_Let1","Pum1_Let2", "Pum2_Let2", "Pum1_Let3")
condensed <- "Comb."
#combine new list into data.frame
pum_let_sum_ex <- data.frame(Group, Fold, condensed, Test)

#make new data frame for grouped pum and let RNA data
pum_let_sum_test <- pum_let_sum
#label values as observed
pum_let_sum_test$Test <- "Observed"
#bind observed and predicted data.frames
pum_let_sum_test <- rbind(pum_let_sum_ex, pum_let_sum_test)

#plot observed and predicted values for let and pum reporters
group1 <- ggplot(subset(pum_let_sum_test, Group == "Pum1_Let1"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-4,0.5)) +
  geom_point(size=2) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank(), legend.position="none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), colour="#56B4E9")
group2 <- ggplot(subset(pum_let_sum_test, Group == "Pum2_Let1"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-4,0.5)) +
  geom_point(size=2) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank(), legend.position="none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), colour="#56B4E9")
group3 <- ggplot(subset(pum_let_sum_test, Group == "Pum3_Let1"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-4,0.5)) +
  geom_point(size=2)  + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank(), legend.position="none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), colour="#56B4E9")
group4 <- ggplot(subset(pum_let_sum_test, Group == "Pum1_Let2"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-4,0.5)) +
  geom_point(size=2)  + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank(), legend.position="none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), colour="#56B4E9")
group5 <- ggplot(subset(pum_let_sum_test, Group == "Pum2_Let2"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-4,0.5)) +
  geom_point(size=2)  + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank(), legend.position="none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), method="loess", colour="#56B4E9")
group6 <- ggplot(subset(pum_let_sum_test, Group == "Pum1_Let3"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-4,0.5)) +
  geom_point(size=2) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank(), legend.position = "none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), colour="#56B4E9")

g_rna_point_pred <-grid.arrange(arrangeGrob(group1, group2, group3, group4, group5, group6, ncol=3, left=textGrob("Relative Expression", rot=90,gp =gpar(fontsize=12, fontfamily="Arial Black"), vjust=0.5)))
ggsave("Pum_Let_RNA_point_pred.tiff", g_rna_point_pred, height=4, width=5, dpi=300, units="in")

ggplot(subset(pum_let_sum_test, Group == "Pum1_Let1"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-4,0.5)) +
  geom_point(size=2) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), colour="#56B4E9")

ggsave("pum_let_point_pred_1.tiff")


#below is the same analysis for TE
TE_results= read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/PARNA_RNA.txt")

Grouped_barcode_identities = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Grouped_Barcode_identities.txt")

TE_results_pum_let <- dplyr::filter(TE_results, RE_Identity %in% Grouped_barcode_identities$Pum_Let_7)

position_pum_let_TE <- subset(TE_results_pum_let, RE_Identity %in% c("7BBB", "7BBp", "7BpB", "7pBB", "pBBB", "pBB7", "pB7B", "p7BB"))

position_pum_let_TE$RE_Identity <- gsub("7", "L", position_pum_let_TE$RE_Identity)
position_pum_let_TE$RE_Identity <- gsub("B", "*", position_pum_let_TE$RE_Identity)

ggplot(position_pum_let_TE, aes(RE_Identity, FOLD)) + 
  geom_boxplot(outlier.colour="white", notch = TRUE) + scale_y_continuous(limits=c(-2.3,1.5)) +
  ylab("Relative TE (log2)") + xlab("RE Identity") + 
  geom_hline(yintercept = 0, linetype="dashed") + theme_science()
ggsave("Pum_Let_position_TE.tiff", height=3, width=5, units="in")

Pum1_Let1_TE <- subset(TE_results_pum_let, condensed %in% c('7p', 'p7', 'p', '7'))
Pum1_Let1_TE$condensed <- gsub("7p|p7", "Comb.", Pum1_Let1_TE$condensed)
Pum2_Let1_TE <- subset(TE_results_pum_let, condensed %in% c('7pp', 'pp7', 'p7p', 'pp', '7'))
Pum2_Let1_TE$condensed <- gsub("7pp|pp7|p7p", "Comb.", Pum2_Let1_TE$condensed)
Pum3_Let1_TE <- subset(TE_results_pum_let, condensed %in% c('7ppp', 'ppp7', 'p7pp', 'pp7p', 'ppp', '7'))
Pum3_Let1_TE$condensed <- gsub("7ppp|ppp7|p7pp|pp7p", "Comb.", Pum3_Let1_TE$condensed)
Pum1_Let2_TE <- subset(TE_results_pum_let, condensed %in% c('77p', 'p77', '7p7', '77', 'p'))
Pum1_Let2_TE$condensed <- gsub("77p|p77|7p7", "Comb.", Pum1_Let2_TE$condensed)
Pum2_Let2_TE <- subset(TE_results_pum_let, condensed %in% c('77pp', '7p7p', '7pp7', 'p77p', 'p7p7', 'pp77', '77', 'pp'))
Pum2_Let2_TE$condensed <- gsub('77pp|7p7p|7pp7|p77p|p7p7|pp77', "Comb.", Pum2_Let2_TE$condensed)
Pum1_Let3_TE <- subset(TE_results_pum_let, condensed %in% c('777p', 'p777', '7p77', '77p7', '777', 'p'))
Pum1_Let3_TE$condensed <- gsub("777p|p777|77p7|7p77", "Comb.", Pum1_Let3_TE$condensed)

Pum1_Let1_TE$Group <- "Pum1_Let1"
Pum2_Let1_TE$Group <- "Pum2_Let1"
Pum3_Let1_TE$Group <- "Pum3_Let1"
Pum1_Let2_TE$Group <- "Pum1_Let2"
Pum2_Let2_TE$Group <- "Pum2_Let2"
Pum1_Let3_TE$Group <- "Pum1_Let3"


pum_let_TE <- rbind(Pum1_Let1_TE, Pum2_Let1_TE, Pum2_Let2_TE, Pum3_Let1_TE, Pum1_Let2_TE, Pum2_Let2_TE, Pum1_Let3_TE)
pum_let_sum_TE <- data.frame(Group = pum_let_TE$Group, Fold = pum_let_TE$FOLD, RE_Identity = pum_let_TE$RE_Identity, condensed = pum_let_TE$condensed)

pum_let_sum_TE$condensed <- gsub("7", "L", pum_let_sum_TE$condensed)

group1 <- ggplot(subset(pum_let_sum_TE, Group == "Pum1_Let1"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("p", "L", "Comb.")) + scale_y_continuous(limits=c(-4,2)) +
  geom_boxplot(outlier.colour="white", fill="white", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 
group2 <- ggplot(subset(pum_let_sum_TE, Group == "Pum2_Let1"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("pp", "L", "Comb.")) + scale_y_continuous(limits=c(-4,2)) +
  geom_boxplot(outlier.colour="white", fill="white", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 
group3 <- ggplot(subset(pum_let_sum_TE, Group == "Pum3_Let1"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("ppp", "L", "Comb.")) + scale_y_continuous(limits=c(-4,2)) +
  geom_boxplot(outlier.colour="white", fill="white", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 
group4 <- ggplot(subset(pum_let_sum_TE, Group == "Pum1_Let2"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("p", "LL", "Comb.")) + scale_y_continuous(limits=c(-4,2)) +
  geom_boxplot(outlier.colour="white", fill="white", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 
group5 <- ggplot(subset(pum_let_sum_TE, Group == "Pum2_Let2"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("pp", "LL", "Comb.")) + scale_y_continuous(limits=c(-4,2)) +
  geom_boxplot(outlier.colour="white", fill="white", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 
group6 <- ggplot(subset(pum_let_sum_TE, Group == "Pum1_Let3"), aes(x = condensed, y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0), limits=c("p", "LLL", "Comb.")) + scale_y_continuous(limits=c(-4,2)) +
  geom_boxplot(outlier.colour="white", fill="white", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 

g_TE <-grid.arrange(arrangeGrob(group1, group2, group3, group4, group5, group6, ncol=3, left=textGrob("Relative TE", rot=90,gp =gpar(fontsize=12, fontfamily="Arial Black"), vjust=0.5)))
ggsave("Pum_Let_TE.tiff", g_TE, height=4, width=5, dpi=300, units="in")

Pum1_Let1_TE <- subset(TE_results_pum_let, condensed %in% c('7p', 'p7', 'p', '7'))
Pum1_Let1_TE$condensed <- gsub("7p|p7", "Comb.", Pum1_Let1_TE$condensed)
Pum1_Let1_TE <- dplyr::group_by(Pum1_Let1_TE, condensed)
Pum1_Let1_TE <- dplyr::summarise(Pum1_Let1_TE, Fold=median(FOLD), mad=mad(FOLD))
Pum2_Let1_TE <- subset(TE_results_pum_let, condensed %in% c('7pp', 'pp7', 'p7p', 'pp', '7'))
Pum2_Let1_TE$condensed <- gsub("7pp|pp7|p7p", "Comb.", Pum2_Let1_TE$condensed)
Pum2_Let1_TE <- dplyr::group_by(Pum2_Let1_TE, condensed)
Pum2_Let1_TE <- dplyr::summarise(Pum2_Let1_TE, Fold=median(FOLD), mad=mad(FOLD))
Pum3_Let1_TE <- subset(TE_results_pum_let, condensed %in% c('7ppp', 'ppp7', 'p7pp', 'pp7p', 'ppp', '7'))
Pum3_Let1_TE$condensed <- gsub("7ppp|ppp7|p7pp|pp7p", "Comb.", Pum3_Let1_TE$condensed)
Pum3_Let1_TE <- dplyr::group_by(Pum3_Let1_TE, condensed)
Pum3_Let1_TE <- dplyr::summarise(Pum3_Let1_TE, Fold=median(FOLD), mad=mad(FOLD))
Pum1_Let2_TE <- subset(TE_results_pum_let, condensed %in% c('77p', 'p77', '7p7', '77', 'p'))
Pum1_Let2_TE$condensed <- gsub("77p|p77|7p7", "Comb.", Pum1_Let2_TE$condensed)
Pum1_Let2_TE <- dplyr::group_by(Pum1_Let2_TE, condensed)
Pum1_Let2_TE <- dplyr::summarise(Pum1_Let2_TE, Fold=median(FOLD), mad=mad(FOLD))
Pum2_Let2_TE <- subset(TE_results_pum_let, condensed %in% c('77pp', '7p7p', '7pp7', 'p77p', 'p7p7', 'pp77', '77', 'pp'))
Pum2_Let2_TE$condensed <- gsub('77pp|7p7p|7pp7|p77p|p7p7|pp77', "Comb.", Pum2_Let2_TE$condensed)
Pum2_Let2_TE <- dplyr::group_by(Pum2_Let2_TE, condensed)
Pum2_Let2_TE <- dplyr::summarise(Pum2_Let2_TE, Fold=median(FOLD), mad=mad(FOLD))
Pum1_Let3_TE <- subset(TE_results_pum_let, condensed %in% c('777p', 'p777', '7p77', '77p7', '777', 'p'))
Pum1_Let3_TE$condensed <- gsub("777p|p777|77p7|7p77", "Comb.", Pum1_Let3_TE$condensed)
Pum1_Let3_TE <- dplyr::group_by(Pum1_Let3_TE, condensed)
Pum1_Let3_TE <- dplyr::summarise(Pum1_Let3_TE, Fold=median(FOLD), mad=mad(FOLD))

Pum1_Let1_TE$Group <- "Pum1_Let1"
Pum2_Let1_TE$Group <- "Pum2_Let1"
Pum3_Let1_TE$Group <- "Pum3_Let1"
Pum1_Let2_TE$Group <- "Pum1_Let2"
Pum2_Let2_TE$Group <- "Pum2_Let2"
Pum1_Let3_TE$Group <- "Pum1_Let3"


pum_let_TE <- rbind(Pum1_Let1_TE, Pum2_Let1_TE, Pum2_Let2_TE, Pum3_Let1_TE, Pum1_Let2_TE, Pum2_Let2_TE, Pum1_Let3_TE)
pum_let_sum_TE <- data.frame(Group = pum_let_TE$Group, Fold = pum_let_TE$Fold, condensed = pum_let_TE$condensed)

pum_let_sum_TE$condensed <- gsub("7", "L", pum_let_sum$condensed)

Test <- "Predicted"
Pum1_Let1_ex_TE <- Pum1_Let1_TE$Fold[Pum1_Let1_TE$condensed == "7"]+Pum1_Let1_TE$Fold[Pum1_Let1_TE$condensed == "p"]
Pum2_Let1_ex_TE <- Pum2_Let1_TE$Fold[Pum2_Let1_TE$condensed == "7"]+Pum2_Let1_TE$Fold[Pum2_Let1_TE$condensed == "pp"]
Pum3_Let1_ex_TE <- Pum3_Let1_TE$Fold[Pum3_Let1_TE$condensed == "7"]+Pum3_Let1_TE$Fold[Pum3_Let1_TE$condensed == "ppp"]
Pum1_Let2_ex_TE <- Pum1_Let2_TE$Fold[Pum1_Let2_TE$condensed == "77"]+Pum1_Let2_TE$Fold[Pum1_Let2_TE$condensed == "p"]
Pum2_Let2_ex_TE <- Pum2_Let2_TE$Fold[Pum2_Let2_TE$condensed == "77"]+Pum2_Let2_TE$Fold[Pum2_Let2_TE$condensed == "pp"]
Pum1_Let3_ex_TE <- Pum1_Let3_TE$Fold[Pum1_Let3_TE$condensed == "777"]+Pum1_Let3_TE$Fold[Pum1_Let3_TE$condensed == "p"]
Fold <- c(Pum1_Let1_ex_TE, Pum2_Let1_ex_TE, Pum3_Let1_ex_TE, Pum1_Let2_ex_TE, Pum2_Let2_ex_TE, Pum1_Let3_ex_TE)
Group <- c("Pum1_Let1", "Pum2_Let1", "Pum3_Let1","Pum1_Let2", "Pum2_Let2", "Pum1_Let3")
condensed <- "Comb."
pum_let_sum_ex <- data.frame(Group, Fold, condensed, Test)

pum_let_sum_test <- pum_let_sum_TE
pum_let_sum_test$Test <- "Observed"

pum_let_sum_test <- rbind(pum_let_sum_ex, pum_let_sum_test)

group1 <- ggplot(subset(pum_let_sum_TE, Group == "Pum1_Let1"), aes(x = reorder(condensed, -Fold), y= 2^Fold)) + 
  xlab("RE Identity") + ylab("Relative TE") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(0.4,1.2)) +
  geom_point(size=2) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank()) + 
  geom_hline(yintercept = 1, linetype="dashed") + geom_smooth(aes(group=1), colour="#56B4E9")
group2 <- ggplot(subset(pum_let_sum_TE, Group == "Pum2_Let1"), aes(x = reorder(condensed, -Fold), y= 2^Fold)) + 
  xlab("RE Identity") + ylab("Relative TE") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(0.4,1.2)) +
  geom_point(size=2) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 1, linetype="dashed") + geom_smooth(aes(group=1), colour="#56B4E9")
group3 <- ggplot(subset(pum_let_sum_TE, Group == "Pum3_Let1"), aes(x = reorder(condensed, -Fold), y= 2^Fold)) + 
  xlab("RE Identity") + ylab("Relative TE") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(0.4,1.2)) +
  geom_point(size=2)  + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 1, linetype="dashed") + geom_smooth(aes(group=1), colour="#56B4E9")
group4 <- ggplot(subset(pum_let_sum_TE, Group == "Pum1_Let2"), aes(x = reorder(condensed, -Fold), y= 2^Fold)) + 
  xlab("RE Identity") + ylab("Relative TE") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(0.4,1.2)) +
  geom_point(size=2)  + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank()) + 
  geom_hline(yintercept = 1, linetype="dashed") + geom_smooth(aes(group=1), colour="#56B4E9")
group5 <- ggplot(subset(pum_let_sum_TE, Group == "Pum2_Let2"), aes(x = reorder(condensed, -Fold), y= 2^Fold)) + 
  xlab("RE Identity") + ylab("Relative TE") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(0.4,1.2)) +
  geom_point(size=2)  + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 1, linetype="dashed") + geom_smooth(aes(group=1), method="loess", colour="#56B4E9")
group6 <- ggplot(subset(pum_let_sum_TE, Group == "Pum1_Let3"), aes(x = reorder(condensed, -Fold), y= 2^Fold)) + 
  xlab("RE Identity") + ylab("Relative TE") + 
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(0.4,1.2)) +
  geom_point(size=2) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 1, linetype="dashed") + geom_smooth(aes(group=1), colour="#56B4E9")

g_TE_point <-grid.arrange(arrangeGrob(group1, group2, group3, group4, group5, group6, ncol=3, left=textGrob("Relative TE", rot=90,gp =gpar(fontsize=12, fontfamily="Arial Black"), vjust=0.5)))
ggsave("Pum_Let_TE_point_linear.tiff", g_TE_point, height=4, width=5, dpi=300, units="in")

group1 <- ggplot(subset(pum_let_sum_test, Group == "Pum1_Let1"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-1,0.5)) +
  geom_point(size=2) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank(), legend.position="none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), colour="#56B4E9")
group2 <- ggplot(subset(pum_let_sum_test, Group == "Pum2_Let1"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-1,0.5)) +
  geom_point(size=2) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank(), legend.position="none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), colour="#56B4E9")
group3 <- ggplot(subset(pum_let_sum_test, Group == "Pum3_Let1"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-1,0.5)) +
  geom_point(size=2)  + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank(), legend.position="none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), colour="#56B4E9")
group4 <- ggplot(subset(pum_let_sum_test, Group == "Pum1_Let2"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-1,0.5)) +
  geom_point(size=2)  + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank(), legend.position="none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), colour="#56B4E9")
group5 <- ggplot(subset(pum_let_sum_test, Group == "Pum2_Let2"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-1,0.5)) +
  geom_point(size=2)  + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank(), legend.position="none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), method="loess", colour="#56B4E9")
group6 <- ggplot(subset(pum_let_sum_test, Group == "Pum1_Let3"), aes(x = reorder(condensed, -Fold), y= Fold, colour=Test)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_y_continuous(limits=c(-1,0.5)) +
  geom_point(size=2) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank(), legend.position = "none") + 
  geom_hline(yintercept = 0, linetype="dashed") + geom_smooth(aes(group=Test), colour="#56B4E9")

g_TE_point_pred <-grid.arrange(arrangeGrob(group1, group2, group3, group4, group5, group6, ncol=3, left=textGrob("Relative TE", rot=90,gp =gpar(fontsize=12, fontfamily="Arial Black"), vjust=0.5)))
ggsave("Pum_Let_TE_point_linear_pred.tiff", g_TE_point_pred, height=4, width=5, dpi=300, units="in")

