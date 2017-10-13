library(dplyr)
library(ggplot2)
library(extrafont)
font_import()

loadfonts(device="win")
cbPallete <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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


RNA_results= read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/small_RNA.txt")

Grouped_barcode_identities = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/Grouped_Barcode_identities.txt")

RNA_results_HuR_Let <- dplyr::filter(RNA_results, RE_Identity %in% Grouped_barcode_identities$HuR_Let_7)

ggplot(RNA_results_HuR_Let, aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + 
  geom_boxplot(notch = TRUE) + theme_science()


Hur1_Let1 <- subset(RNA_results_HuR_Let, condensed %in% c('7h', 'h7', 'h','7', 'control'))

Hur1_Let1$Group <- ""

Hur1_Let1$Group1[Hur1_Let1$RE_Identity %in%  c("7BBB", "BhBB", "7hBB")] <- 1
Hur1_Let1$Group2[Hur1_Let1$RE_Identity %in%  c("7BBB", "BBhB", "7BhB")] <- 2
Hur1_Let1$Group3[Hur1_Let1$RE_Identity %in%  c("7BBB", "BBBh", "7BBh")] <- 3
Hur1_Let1$Group4[Hur1_Let1$RE_Identity %in%  c("B7BB", "hBBB", "h7BB")] <- 4
Hur1_Let1$Group5[Hur1_Let1$RE_Identity %in%  c("B7BB", "BBhB", "B7hB")] <- 5
Hur1_Let1$Group6[Hur1_Let1$RE_Identity %in%  c("B7BB", "BBBh", "B7Bh")] <- 6
Hur1_Let1$Group7[Hur1_Let1$RE_Identity %in%  c("BB7B", "hBBB", "hB7B")] <- 7
Hur1_Let1$Group8[Hur1_Let1$RE_Identity %in%  c("BB7B", "BhBB", "Bh7B")] <- 8
Hur1_Let1$Group9[Hur1_Let1$RE_Identity %in%  c("BB7B", "BBBh", "BB7h")] <- 9
Hur1_Let1$Group10[Hur1_Let1$RE_Identity %in%  c("BBB7", "hBBB", "hBB7")] <- 10
Hur1_Let1$Group11[Hur1_Let1$RE_Identity %in%  c("BBB7", "BhBB", "BhB7")] <- 11
Hur1_Let1$Group12[Hur1_Let1$RE_Identity %in%  c("BBB7", "BBhB", "BBh7")] <- 12

Hur1_Let1$RE_Identity <- gsub("B", "*", Hur1_Let1$RE_Identity)
Hur1_Let1$RE_Identity <- gsub("7", "L", Hur1_Let1$RE_Identity)
Hur1_Let1$RE_Identity <- gsub("h", "A", Hur1_Let1$RE_Identity)

group1 <- ggplot(subset(Hur1_Let1, Group1 == 1), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour="white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed") 
group2 <- ggplot(subset(Hur1_Let1, Group2 == 2), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour="white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed")
group3 <- ggplot(subset(Hur1_Let1, Group3 == 3), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour="white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank())+ 
  geom_hline(yintercept = 0, linetype="dashed")
group4 <- ggplot(subset(Hur1_Let1, Group4 == 4), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour="white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed")
group5 <- ggplot(subset(Hur1_Let1, Group5 == 5), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour="white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed")
group6 <- ggplot(subset(Hur1_Let1, Group6 == 6), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour="white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed")
group7 <- ggplot(subset(Hur1_Let1, Group7 == 7), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour="white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed")
group8 <- ggplot(subset(Hur1_Let1, Group8 == 8), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour="white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed")
group9 <- ggplot(subset(Hur1_Let1, Group9 == 9), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour="white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed")
group10 <- ggplot(subset(Hur1_Let1, Group10 == 10), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour="white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed")
group11 <- ggplot(subset(Hur1_Let1, Group11 == 11), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour="white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed")
group12 <- ggplot(subset(Hur1_Let1, Group12 == 12), aes(x = RE_Identity, y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + scale_y_continuous(limits=c(-3.1,3.1)) +
  geom_boxplot(outlier.colour = "white", fill="#009E73", notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  geom_hline(yintercept = 0, linetype="dashed")

g_rna <-grid.arrange(arrangeGrob(group1, group2, group3, group4, group5, group6, group7, group8, group9, group10, group11, group12, ncol=6, left=textGrob("Relative 40s Association (log2)", rot=90,gp =gpar(fontsize=12, fontfamily="Arial Black"), vjust=0.5)))
ggsave("Hur1_Let1_small.tiff", g_rna, height=4, width=5, dpi=300, units="in")

