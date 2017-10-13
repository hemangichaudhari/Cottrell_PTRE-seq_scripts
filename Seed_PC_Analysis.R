library(dplyr)
library(ggplot2)
library(reshape2)
library(extrafont)
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

bar_order <- c("Control", "6mer_x1", "6mer_x2", '6mer_x4', "7mer_a1_x1", "7mer_a1_x2", "7mer_a1_x4", "7mer_m8_x1", "7mer_m8_x2","7mer_m8_x4", "8mer_x1", "8mer_x2", "8mer_x4", "PC_x1", "PC_x2", "PC_x4")

small_results= read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/small_RNA.txt")

small_results_seed <- small_results

Grouped_barcode_identities = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/Grouped_Barcode_identities.txt")

small_results_seed <- dplyr::filter(small_results_seed, RE_Identity %in% Grouped_barcode_identities$Let_PC)

small_results_seed$RE_Identity <- gsub("7777", "8mer_x4", small_results_seed$RE_Identity )
small_results_seed$RE_Identity <- gsub("7B7B", "8mer_x2", small_results_seed$RE_Identity )
small_results_seed$RE_Identity <- gsub("B7BB", "8mer_x1",small_results_seed$RE_Identity )

small_results_seed$RE_Identity <- gsub("7mera1x", "7mer_a1_x", small_results_seed$RE_Identity)

small_results_seed$RE_Identity <- gsub("7merm8x", "7mer_m8_x", small_results_seed$RE_Identity)

small_results_seed$RE_Identity <- gsub("6merx", "6mer_x", small_results_seed$RE_Identity)

small_results_seed$RE_Identity <- gsub("7pcx", "PC_x", small_results_seed$RE_Identity)

bar_order <- c("Control", "6mer_x1", "6mer_x2", '6mer_x4', "7mer_a1_x1", "7mer_a1_x2", "7mer_a1_x4", "7mer_m8_x1", "7mer_m8_x2","7mer_m8_x4", "8mer_x1", "8mer_x2", "8mer_x4", "PC_x1", "PC_x2", "PC_x4")

ggplot(small_results_seed, aes(RE_Identity, FOLD)) + 
  xlab("") + ylab("Relative 40s Assoc.(log2)") +  scale_x_discrete(limits=bar_order) +
  geom_boxplot(outlier.shape = 1) + theme_science() + geom_hline(yintercept = 0, linetype="dashed") + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank())

ggsave("Seed_small.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

seed_pc <- ggplot(small_results_seed, aes(RE_Identity, FOLD)) + 
  xlab("") + ylab("Relative 40s Assoc.(log2)") +  scale_x_discrete(limits=bar_order) +
  geom_boxplot(outlier.colour = NA) + theme_science() +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")

sts_pc_1 <- boxplot.stats(small_results_seed$FOLD)$stats
p1 = seed_pc + coord_cartesian(ylim = c(min(sts_pc_1), max(sts_pc_1)))
p1

ggsave("Seed_small_NO.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#t.test(small_results$FOLD[small_results$condensed=="p"], small_results$FOLD[small_results$condensed=="Control"])
#t.test(small_results$FOLD[small_results$condensed=="pp"], small_results$FOLD[small_results$condensed=="Control"])
#t.test(small_results$FOLD[small_results$condensed=="ppp"], small_results$FOLD[small_results$condensed=="Control"])
#t.test(small_results$FOLD[small_results$condensed=="pppp"], small_results$FOLD[small_results$condensed=="Control"])

monosome_results= read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/monosome_RNA.txt")

monosome_results_seed <- monosome_results

Grouped_barcode_identities = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/Grouped_Barcode_identities.txt")

monosome_results_seed <- dplyr::filter(monosome_results_seed, RE_Identity %in% Grouped_barcode_identities$Let_PC)

monosome_results_seed$RE_Identity <- gsub("7777", "8mer_x4", monosome_results_seed$RE_Identity )
monosome_results_seed$RE_Identity <- gsub("7B7B", "8mer_x2", monosome_results_seed$RE_Identity )
monosome_results_seed$RE_Identity <- gsub("B7BB", "8mer_x1",monosome_results_seed$RE_Identity )

monosome_results_seed$RE_Identity <- gsub("7mera1x", "7mer_a1_x", monosome_results_seed$RE_Identity)

monosome_results_seed$RE_Identity <- gsub("7merm8x", "7mer_m8_x", monosome_results_seed$RE_Identity)

monosome_results_seed$RE_Identity <- gsub("6merx", "6mer_x", monosome_results_seed$RE_Identity)

monosome_results_seed$RE_Identity <- gsub("7pcx", "PC_x", monosome_results_seed$RE_Identity)


ggplot(monosome_results_seed, aes(RE_Identity, FOLD)) + 
  xlab("") + ylab("Relative MA (log2)") +  scale_x_discrete(limits=bar_order) +
  geom_boxplot(outlier.shape = 1) + theme_science() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")

ggsave("Seed_monosome.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

seed_pc_te <- ggplot(monosome_results_seed, aes(RE_Identity, FOLD)) + 
  xlab("") + ylab("Relative MA (log2)") +  scale_x_discrete(limits=bar_order) +
  geom_boxplot(outlier.colour = NA) + theme_science() +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")

sts_pc_2 <- boxplot.stats(small_results_seed$FOLD)$stats
p2 = seed_pc_te + coord_cartesian(ylim = c(min(sts_pc_2), max(sts_pc_2)))
p2

ggsave("Seed_monosome_NO.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

