library(dplyr)
library(ggplot2)
library(extrafont)
library(reshape2)
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

#make list for bar_order
bar_order = c("Control", "1", "2", "3", "4")
#load RNA results
RNA_results= read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/RNA_plasmid.txt")
#subset RNA results for only Let-7 targeted reporters
RNA_results_let_7 <- subset(RNA_results, condensed %in% c("7", "77", "777", "7777", "Control"))
#rename condensed by number of Let-7 sites
RNA_results_let_7$condensed <- gsub("7777", "4", RNA_results_let_7$condensed)                            
RNA_results_let_7$condensed <- gsub("777", "3", RNA_results_let_7$condensed)                                
RNA_results_let_7$condensed <- gsub("77", "2", RNA_results_let_7$condensed)    
RNA_results_let_7$condensed <- gsub("7", "1", RNA_results_let_7$condensed) 
#rename B -> * and 7 -> L
RNA_results_let_7$RE_Identity <- gsub("B", "*", RNA_results_let_7$RE_Identity)  
RNA_results_let_7$RE_Identity <- gsub("7", "L", RNA_results_let_7$RE_Identity)  
#plot RNA expression by number of Let-7 sites
ggplot(RNA_results_let_7, aes(condensed, FOLD, fill=condensed)) + scale_fill_manual(guide=FALSE, values = cbPalette) + 
  xlab("Number of Let-7 sites") + ylab("Relative Expression (log2)") +  scale_x_discrete(limits=bar_order) +
  geom_boxplot(outlier.shape = 1) + theme_science() + geom_hline(yintercept = 0, linetype="dashed")

ggsave("Let-7_RNA.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

ggplot(RNA_results_let_7, aes(condensed, FOLD, fill=condensed)) + scale_fill_manual(guide=FALSE, values = cbPalette) + 
  xlab("Number of Let-7 sites") + ylab("Relative Expression (log2)") +  scale_x_discrete(limits=bar_order) +
  geom_jitter() + theme_science() + geom_hline(yintercept = 0, linetype="dashed")

#plot RNA expression by number of Let-7 sites with invisible outliers
let_7 <- ggplot(RNA_results_let_7, aes(condensed, FOLD, fill=condensed)) + scale_fill_manual(guide=FALSE, values = cbPalette) + 
  xlab("Number of Let-7 sites") + ylab("Relative Expression (log2)") +  scale_x_discrete(limits=bar_order) +
  geom_boxplot(outlier.colour =  NA, notch = TRUE) + theme_science() + geom_hline(yintercept = 0, linetype="dashed")

#rescale plot form above
sts_1 <- boxplot.stats(RNA_results_let_7$FOLD)$stats
p1 = let_7 + coord_cartesian(ylim = c(min(sts_1),max(sts_1)))
p1

ggsave("Let-7_RNA_NO.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#remake above plot colored orange
let_7 <- ggplot(RNA_results_let_7, aes(condensed, FOLD)) + 
  xlab("Number of Let-7 sites") + ylab("Relative Expression (log2)") +  scale_x_discrete(limits=bar_order) +
  annotate("text", x=c("Control", "1", "2", "3", "4"), y=c(1.7,0.5,-0.4,-0.35,-1.1), fontface=2, hjust=0.4, size=4, label=(c("", "***", "***", "***", "***"))) +
  geom_boxplot(outlier.colour =  NA, fill="#E69F00", notch = TRUE) + theme_science() + geom_hline(yintercept = 0, linetype="dashed") 


sts_1 <- boxplot.stats(RNA_results_let_7$FOLD)$stats
p1 = let_7 +  coord_cartesian(ylim = c(min= -4.2,max=2.6)) + labs(title="Let-7") + theme(plot.title = element_text(hjust = 0.5, size=24))
p1

ggsave("Let-7_RNA_NO_orange.tiff", plot=last_plot(), width=5, height=3.5, units="in", dpi=300)

#plot RNA expression of let-7 reporters
ggplot(RNA_results_let_7, aes(x = reorder(RE_Identity, -FOLD), FOLD, fill=condensed)) + scale_fill_manual(name = "Number of Sites", values = cbPalette) +
  xlab("") + ylab("Relative Expression (log2)") + 
  geom_boxplot(outlier.shape = 1, notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")

ggsave("Let-7_RNA_All.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

ggplot(RNA_results_let_7, aes(x = reorder(RE_Identity, -FOLD), FOLD, fill=condensed)) + scale_colour_manual(name = "Number of Sites", values = cbPalette) +
  xlab("") + ylab("Relative Expression (log2)") + 
  geom_jitter() + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")

#plot RNA expression of let-7 reporters with invisible outliers
let_7_all <- ggplot(RNA_results_let_7, aes(x = reorder(RE_Identity, -FOLD), FOLD, fill=condensed)) + scale_fill_manual(name = "Number of Sites", values = cbPalette) +
  xlab("") + ylab("Relative Expression (log2)") +  
  geom_boxplot(outlier.colour = NA, notch = TRUE) + theme_science() + theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")
#rescale plot from above
sts_2 <- boxplot.stats(RNA_results_let_7$FOLD)$stats
p2 = let_7_all + coord_cartesian(ylim = c(min(sts_2),max(sts_2)))
p2

ggsave("Let-7_RNA_All_NO.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#load TE results
TE_results= read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/PARNA_RNA.txt")
#subset TE results for only let-7 targeted reporters
TE_results_let_7 <- subset(TE_results, condensed %in% c("7", "77", "777", "7777", "Control"))
#rename condensed by number of let-7 sites
TE_results_let_7$condensed <- gsub("7777", "4", TE_results_let_7$condensed)                            
TE_results_let_7$condensed <- gsub("777", "3", TE_results_let_7$condensed)                                
TE_results_let_7$condensed <- gsub("77", "2", TE_results_let_7$condensed)    
TE_results_let_7$condensed <- gsub("7", "1", TE_results_let_7$condensed)  
#remane B -> * and 7 -> L
TE_results_let_7$RE_Identity <- gsub("B", "*", TE_results_let_7$RE_Identity)  
TE_results_let_7$RE_Identity <- gsub("7", "L", TE_results_let_7$RE_Identity)  
#plot TE and number of let-7 sites
ggplot(TE_results_let_7, aes(condensed, FOLD, fill=condensed)) + scale_fill_manual(guide=FALSE, values = cbPalette) +
  xlab("Number of Let-7 sites") + ylab("Relative TE (log2)") +  scale_x_discrete(limits=bar_order) +
  geom_boxplot(outlier.shape = 1, notch = TRUE) + theme_science() + geom_hline(yintercept = 0, linetype="dashed")

ggsave("Let-7_TE.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)
#plot TE and number of Let-7 sites with invisible outliers
let_7_TE <- ggplot(TE_results_let_7, aes(condensed, FOLD, fill=condensed)) + scale_fill_manual(guide=FALSE, values = cbPalette) +
  xlab("Number of Let-7 sites") + ylab("Relative TE (log2)") +  scale_x_discrete(limits=bar_order) +
  geom_boxplot(outlier.colour = NA, notch = TRUE) + theme_science() + geom_hline(yintercept = 0, linetype="dashed")
#rescale above plot
sts_3 <- boxplot.stats(TE_results_let_7$FOLD)$stats
p3 = let_7_TE + coord_cartesian(ylim = c(min(sts_3)*1.5,max(sts_3)))
p3

ggsave("Let-7_TE_NO.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)
#remake above plot colored blue
let_7_TE <- ggplot(TE_results_let_7, aes(condensed, FOLD)) +
  xlab("Number of Let-7 sites") + ylab("Relative TE (log2)") +  scale_x_discrete(limits=bar_order) +
  annotate("text", x=c("Control", "1", "2", "3", "4"), y=c(1.7,1.1,1.3,1.1,1), fontface=2, hjust=0.4, size=4, label=(c("", "***", "***", "***", "***"))) +
  geom_boxplot(outlier.colour = NA, fill="#56B4E9", notch = TRUE) + theme_science() + geom_hline(yintercept = 0, linetype="dashed")

sts_3 <- boxplot.stats(TE_results_let_7$FOLD)$stats
p3 = let_7_TE + coord_cartesian(ylim = c(min= -4.2,max=2.6))
p3

ggsave("Let-7_TE_NO_blue.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)
#plot TE for all let-7 targeted reporters
ggplot(TE_results_let_7, aes(x = reorder(RE_Identity, -FOLD), FOLD, fill=condensed)) + scale_fill_manual(name = "Number of Sites", values = cbPalette) +
  xlab("") + ylab("Relative TE (log2)") + 
  geom_boxplot(outlier.shape = 1, notch = TRUE) + theme_science() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")

ggsave("Let-7_TE_All.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)
#plot TE for all let-7 targeted reporters, jitter
ggplot(TE_results_let_7, aes(x = reorder(RE_Identity, -FOLD), FOLD, fill=condensed)) + scale_fill_manual(name = "Number of Sites", values = cbPalette) +
  xlab("") + ylab("Relative TE (log2)") + 
  geom_jitter() + theme_science() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")

#plot TE for all let-7 reporters with outliers invisible
let_7_TE_all <- ggplot(TE_results_let_7, aes(x = reorder(RE_Identity, -FOLD), FOLD, fill=condensed)) + scale_fill_manual(name = "Number of Sites", values = cbPalette) +
  xlab("") + ylab("Relative TE (log2)") + 
  geom_boxplot(outlier.colour=NA, notch = TRUE) + theme_science() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept = 0, linetype="dashed")
#rescale above plot
sts_4 <- boxplot.stats(TE_results_let_7$FOLD)$stats
p4 = let_7_TE_all + coord_cartesian(ylim = c(min(sts_4)*1.5,max(sts_4)))
p4

ggsave("Let-7_TE_All_NO.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#perform statistics for let-7 effects on RNA
summary(aov(FOLD ~ condensed, data=RNA_results_let_7))
let_rna_stats <- pairwise.t.test(RNA_results_let_7$FOLD, RNA_results_let_7$condensed, p.adj = "bonf")
let_rna_stats <- let_rna_stats$p.value
let_rna_stats[let_rna_stats < 2e-16 ] <- "< 2e-16"
let_rna_stats <- melt(let_rna_stats)
let_rna_stats <- na.omit(let_rna_stats)
colnames(let_rna_stats) <- c("Var1", "Var2", "RNA.p.value")


#perform statistics for let-7 effects on TE
summary(aov(FOLD ~ condensed, data=TE_results_let_7))
let_TE_stats <- pairwise.t.test(TE_results_let_7$FOLD, TE_results_let_7$condensed, p.adj = "bonf")
let_TE_stats <- let_TE_stats$p.value
let_TE_stats[let_TE_stats < 2e-16 ] <- "< 2e-16"
let_TE_stats <- melt(let_TE_stats)
let_TE_stats <- na.omit(let_TE_stats)
colnames(let_TE_stats) <- c("Var1", "Var2", "TE.p.value")

#merge t-trest results and write table
let_stats <- merge(let_TE_stats, let_rna_stats)

write.csv(let_stats, "Let_t_tests.csv")

#assign new column as level (TE or RNA)
TE_results_let_7$Level <- "TE"
RNA_results_let_7$Level <- "RNA"

#bind TE and RNA results
results_let_7 <- rbind(TE_results_let_7, RNA_results_let_7)

#plot TE and RNA by number of let-7 binding sites
ggplot(results_let_7, aes(condensed, FOLD, fill=Level)) + scale_fill_manual(values = c("gray", "white")) +
  xlab("Number of Let-7 sites") + ylab("Fold Change (log2)") +  scale_x_discrete(limits=bar_order) +
  geom_boxplot(outlier.shape = 1, notch = TRUE) + theme_science() + geom_hline(yintercept = 0, linetype="dashed")

ggsave("Let-7_TE_RNA.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#plot TE and RNA by number of let-7 binding sites with invisible outliers
let_7_TE_RNA <- ggplot(results_let_7, aes(condensed, FOLD, fill=Level)) + scale_fill_manual(values = c("gray", "white")) +
  xlab("Number of Let-7 sites") + ylab("Fold Change (log2)") +  scale_x_discrete(limits=bar_order) +
  geom_boxplot(outlier.colour = NA, notch = TRUE) + theme_science() + geom_hline(yintercept = 0, linetype="dashed") +
  annotate("text", x=c("Control", "1", "2", "3", "4"), y=1.3, fontface=4, hjust=-0.15, size=4, label=(c("", "***", "***", "***", "***"))) +
  annotate("text", x=c("Control", "1", "2", "3", "4"), y=c(1, 0.5,-0.25, -0.2, -1), colour="grey", fontface=4, hjust=1.2, size=4, label=(c("", "***", "***", "***", "***")))
  
#rescale above plot
sts_5 <- boxplot.stats(results_let_7$FOLD)$stats
p5 = let_7_TE_RNA + coord_cartesian(ylim = c(min(sts_5)*1.1,max(sts_5)))
p5


ggsave("Let-7_TE_RNA_NO.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)

#plot TE and RNA by number of let-7 binding sites on a linear scale
let_7_TE_RNA_linear <- ggplot(results_let_7, aes(condensed, 2^FOLD, fill=Level)) + scale_fill_manual(values = c("gray", "white")) +
  xlab("Number of Let-7 sites") + ylab("Fold Change") +  scale_x_discrete(limits=bar_order) +
  geom_smooth(aes(group=Level, colour=Level), se=FALSE) + scale_color_manual(values=c("grey", "black")) +
  geom_boxplot(outlier.colour = NA, notch = TRUE) + theme_science() + geom_hline(yintercept = 1, linetype="dashed") 

sts_6 <- boxplot.stats(2^results_let_7$FOLD)$stats
p6 = let_7_TE_RNA_linear + coord_cartesian(ylim = c(min(sts_6),max(sts_6)))
p6

ggsave("Let-7_TE_RNA_NO_linear.tiff", plot=last_plot(), width=5, height=3, units="in", dpi=300)
