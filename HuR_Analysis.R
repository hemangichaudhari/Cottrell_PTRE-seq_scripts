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

small_results= read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/small_RNA.txt")
small_results_HuR <- subset(small_results, condensed %in% c("h", "hh", "hhh", "hhhh", "Control"))

small_results_HuR$condensed <- gsub("hhhh", "4", small_results_HuR$condensed)                            
small_results_HuR$condensed <- gsub("hhh", "3", small_results_HuR$condensed)                                
small_results_HuR$condensed <- gsub("hh", "2", small_results_HuR$condensed)    
small_results_HuR$condensed <- gsub("h", "1", small_results_HuR$condensed)  

small_results_HuR$RE_Identity <- gsub("B", "*", small_results_HuR$RE_Identity)
small_results_HuR$RE_Identity <- gsub("h", "A", small_results_HuR$RE_Identity)

ggplot(small_results_HuR, aes(x = RE_Identity, FOLD, fill=condensed)) + scale_fill_manual(name = "Number of Sites", values = cbPalette) + 
  xlab("RE Identity") + ylab("Relative 40s Assoc. (log2)") + geom_hline(yintercept = 0, linetype="dashed") +
  geom_boxplot(outlier.shape = 1, notch = TRUE) + theme_science() + scale_x_discrete(limits=c("AAA*", "AA**", "A*A*", "*AA*", "*A**", "Control", "**A*", "A***", "***A", "A**A", "AA*A","*A*A", "A*AA", "**AA", "AAAA", "*AAA")) +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank())

ggsave("HuR_small.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

HuR_1 <- ggplot(small_results_HuR, aes(x = reorder(RE_Identity, FOLD), FOLD)) + 
  xlab("RE Identity") + ylab("Relative 40s Assoc. (log2)") + geom_hline(yintercept = 0, linetype="dashed") +
  geom_boxplot(outlier.colour = NA, fill="#009E73", notch = TRUE) + theme_science() + scale_x_discrete(limits=c("AAA*", "AA**", "A*A*", "*AA*", "*A**", "Control", "**A*", "A***", "***A", "A**A", "AA*A","*A*A", "A*AA", "**AA", "AAAA", "*AAA")) +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank())

sts_HuR_1 <- boxplot.stats(small_results_HuR$FOLD)$stats
p1 = HuR_1 + coord_cartesian(ylim = c(max(sts_HuR_1)*0.95,min(sts_HuR_1)*0.7))
p1

ggsave("HuR_small_NO.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)


monosome_results= read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_Ribosome/monosome_RNA.txt")

monosome_results_HuR <- subset(monosome_results, condensed %in% c("h", "hh", "hhh", "hhhh", "Control"))

monosome_results_HuR$condensed <- gsub("hhhh", "4", monosome_results_HuR$condensed)                            
monosome_results_HuR$condensed <- gsub("hhh", "3", monosome_results_HuR$condensed)                                
monosome_results_HuR$condensed <- gsub("hh", "2", monosome_results_HuR$condensed)    
monosome_results_HuR$condensed <- gsub("h", "1", monosome_results_HuR$condensed)  



ggplot(monosome_results_HuR, aes(x = reorder(RE_Identity, -FOLD), FOLD, fill=condensed)) + scale_fill_manual(name = "Number of Sites", values = cbPalette) +
  xlab("RE Identity") + ylab("Relative MA (log2)") + geom_hline(yintercept = 0, linetype="dashed") +
  geom_boxplot(outlier.shape = 1) + theme_science() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank())

ggsave("HuR_monosome.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

HuR_2 <- ggplot(monosome_results_HuR, aes(x = reorder(RE_Identity, -FOLD), FOLD, fill=condensed)) + scale_fill_manual(name = "Number of Sites", values = cbPalette) + 
  xlab("RE Identity") + ylab("Relative MA (log2)") + geom_hline(yintercept = 0, linetype="dashed") +
  geom_boxplot(outlier.colour = NA) + theme_science() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank())

sts_HuR_2 <- boxplot.stats(monosome_results_HuR$FOLD)$stats
p2 = HuR_2 + coord_cartesian(ylim = c(max(sts_HuR_2)*0.8,min(sts_HuR_2)*1.5))
p2

ggsave("HuR_monosome_NO.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

