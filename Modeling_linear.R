# Kyle's microRNA and RBP project
#3 Feb 2017
#Hemangi Chaudhari

#---Call libraries----------------------------------------------------------------------------------------
library(seqinr)
library(ggplot2)
library(reshape2)
library(pROC)
library(seqinr)
library(ggplot2)
library("gkmSVM")
library("gtools")
library("boot")
library(leaps)
library(dplyr)

library(extrafont)
library(gridExtra)
library(RGraphics)
library(pwr)

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



#---Set Working Directory (Input Output Files)------------------------------------------------------------
setwd("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1")

#---INPUT-------------------------------------------------------------------------------------------------
bcid = read.delim("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/Grouped_Barcode_identities.txt")
TE = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/PARNA_RNA.txt",header = TRUE)
RNA = read.table("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/RNA_plasmid.txt",header = TRUE)
Seq = read.table("/Users/Kyle/Dropbox/CRE_seq/Library.txt",header = TRUE)

Seq$RE_Identity <- gsub("-", "B", Seq$RE_Identity)

#---Re-format Inputs------------------------------------------------------------------------------------
# Keep unique sequences without barcodes
Seq = Seq[,-c(2)] # Removes second column
Seq = unique(Seq) # Keeps only unique rows
# Calculate mean FOLD expression for all barcodes for a given pattern
aggregate(FOLD ~ RE_Identity, data= RNA, FUN=median) -> RNA_median
aggregate(FOLD ~ RE_Identity, data= RNA, FUN=mean) -> RNA_mean
# Add Sequence information to pattern  
RNA_median = merge(RNA_median, Seq, by ="RE_Identity", all.x = TRUE)
RNA = RNA_median
# Change control to BBBB
RNA[,1] <- as.character(RNA[,1])
RNA[which(RNA[,1] == "Control"),1] <- ("BBBB")

colnames(RNA) <- c("RE_Identity","FOLD","sequence")

# Split the pattern to get individual positions 
X =(strsplit(as.character(RNA$RE_Identity), split= ""))
n.obs <- sapply(X, length) # deals with different column numbers in each row
seq.max <- seq_len(max(n.obs))# deals with different column numbers in each row
X <- t(sapply(X, "[", i = seq.max))# deals with different column numbers in each row
X[X == 7] <- "L" # Replace 7 with L 
colnames(X) = paste("P", c(1:8), sep ="") #Change column names 
RNA = cbind(RNA, X) 
RNA$type = nchar((as.character(RNA$RE_Identity)))
RNA_4 = subset(RNA, type == 4) # Keep guys will 4 position pattern only
RNA_4 = subset(RNA_4, RE_Identity != "Dna2" ) # Remove Dna2 guy
#RNA_4<- subset(RNA_4, !grepl("h", RNA_4$RE_Identity))

#---Linear Models-----------------------------------------------------------------------------------
# Set Reference level to B for linear model fits 
RNA_4 <- within(RNA_4, P1 <- relevel(P1, ref = "B"))
RNA_4 <- within(RNA_4, P2 <- relevel(P2, ref = "B"))
RNA_4 <- within(RNA_4, P3 <- relevel(P3, ref = "B"))
RNA_4 <- within(RNA_4, P4 <- relevel(P4, ref = "B"))

# Formula for the linear model 
fmla <- as.formula(paste("FOLD ~ P1*P2 + P1*P3 + P1*P4 + P2*P3 + P2*P4 + P3*P4  "))
# this particular model tries to predict expression with all 4 positions and all interactions
# P1*P2 + P1*P3 + P1*P4 + P2*P3 + P2*P4 + P3*P4
# Fit all data
fit <-  glm(fmla, data = RNA_4)
summary(fit) #displays summary to the fit. Information stored in this variable
predpr<- predict(fit,type=c("response")) # Predicted values from the fit 
RNA_4$problm=predpr # Adding predicted values to the dataset
t= cor(RNA_4$problm, RNA_4$FOLD)
ggplot(RNA_4, aes(FOLD,problm)) + geom_point() + 
  ggtitle(label=(paste("R =", format(t, digits = 2)))) + 
  xlab("Measured Relative RNA Expression (log2)") + 
  ylab("Predicted Relative RNA Expression (log2)") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()
ggsave("Fit_all_interactions_RNA.tiff", width=5, height = 5, unit = "in")

f2_rna= (t^2)/(1-t^2)

pwr.f2.test(u=112, v=512, f2= f2_rna, sig.level = 0.01, power=NULL)


#model predictions for onyl let-7 targets filter RNA_4 and find only let-7 targets
let_7_RNA <- dplyr::filter(RNA_4, RE_Identity %in% bcid_rename$Let.7)
#correlation between predicted and observed
t_let= cor.test(let_7_RNA$problm, let_7_RNA$FOLD)
#plot predicted v observed
ggplot(let_7_RNA, aes(FOLD,problm)) + geom_point() + 
  ggtitle (paste("R=", round(t_let$estimate, digits = 3))) + 
  xlab("Measured Relative RNA Expression (log2)") + 
  ylab("Predicted Relative RNA Expression (log2)") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()
ggsave("model_prediction_let_rna.tiff", height=6, width=6, units="in")

#model predictions for onyl let-7 targets filter RNA_4 and find only Pumilio targets
pum_RNA <- dplyr::filter(RNA_4, RE_Identity %in% bcid_rename$Pumilio)
#correlation between predicted and observed
t_pum= cor.test(pum_RNA$problm, pum_RNA$FOLD)
#plot predicted v observed
ggplot(pum_RNA, aes(FOLD,problm)) + geom_point() + 
  ggtitle (paste(i,": R=", round(t_pum$estimate, digits = 3))) + 
  xlab("Measured Relative RNA Expression") + 
  ylab("Predicted Relative RNA Expression") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()

#model predictions for onyl let-7 targets filter RNA_4 and find only HuR (ARE) targets
hur_RNA <- dplyr::filter(RNA_4, RE_Identity %in% bcid_rename$HuR)
#correlation between predicted and observed
t_hur= cor.test(hur_RNA$problm, hur_RNA$FOLD)
#plot predicted v observed
ggplot(let_7_RNA, aes(FOLD,problm)) + geom_point() + 
  ggtitle (paste(i,": R=", round(t_hur$estimate, digits = 3))) + 
  xlab("Measured Relative RNA Expression") + 
  ylab("Predicted Relative RNA Expression") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()

#model predictions for onyl let-7 targets filter RNA_4 and find only Smg targets
smg_RNA <- dplyr::filter(RNA_4, RE_Identity %in% bcid_rename$Smg)
#correlation between predicted and observed
t_smg= cor.test(smg_RNA$problm, smg_RNA$FOLD)
#plot predicted v observed
ggplot(smg_RNA, aes(FOLD,problm)) + geom_point() + 
  ggtitle (paste(i,": R=", round(t_smg$estimate, digits = 3))) + 
  xlab("Measured Relative RNA Expression") + 
  ylab("Predicted Relative RNA Expression") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()

#model predictions for onyl let-7 targets filter RNA_4 and find only Pum+Let targets
pum_let_RNA <- dplyr::filter(RNA_4, RE_Identity %in% bcid_rename$Pum_Let_7)
#correlation between predicted and observed
t_pum_let= cor.test(pum_let_RNA$problm, pum_let_RNA$FOLD)
#plot predicted v observed
ggplot(pum_let_RNA, aes(FOLD,problm)) + geom_point() + 
  ggtitle (paste(i,": R=", round(t_pum_let$estimate, digits = 3))) + 
  xlab("Measured Relative RNA Expression") + 
  ylab("Predicted Relative RNA Expression") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()

#model predictions for onyl let-7 targets filter RNA_4 and find only Pum+ARE targets
pum_hur_RNA <- dplyr::filter(RNA_4, RE_Identity %in% bcid_rename$HuR_Pum)
#correlation between predicted and observed
t_pum_hur= cor.test(pum_hur_RNA$problm, pum_hur_RNA$FOLD)
#plot predicted v observed
ggplot(pum_hur_RNA, aes(FOLD,problm)) + geom_point() + 
  ggtitle (paste(i,": R=", round(t_pum_hur$estimate, digits = 3))) + 
  xlab("Measured Relative RNA Expression") + 
  ylab("Predicted Relative RNA Expression") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()

#model predictions for onyl let-7 targets filter RNA_4 and find only Let-7+ARE targets
let_hur_RNA <- dplyr::filter(RNA_4, RE_Identity %in% bcid_rename$HuR_Let_7)
#correlation between predicted and observed
t_let_hur= cor.test(let_hur_RNA$problm, let_hur_RNA$FOLD)
#plot predicted v observed
ggplot(let_hur_RNA, aes(FOLD,problm)) + geom_point() + 
  ggtitle (paste(i,": R=", round(t_let_hur$estimate, digits = 3))) + 
  xlab("Measured Relative RNA Expression") + 
  ylab("Predicted Relative RNA Expression") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()


data.frame(summary(fit)$coeff)-> coeff # stores information about paramters, estimates and significance 
coeff_RNA=coeff[with(coeff, order(coeff[,4])), ] # you can order coeffs by Effect (Column1) or Signficance (Column 4)


# Perform fit 5 fold cross validation
RNA_4<-RNA_4[sample(nrow(RNA_4)),]  # Randomize rows
folds <- cut(seq(1,nrow(RNA_4)),breaks=5,labels=FALSE) #Assigns every row to be one part of 5 folds

for(i in 1:5){
  testIndexes <- which(folds==i,arr.ind=TRUE) #Segement your data by fold using the which() function 
  testData <-RNA_4[testIndexes, ] # Data assinged to current index
  trainData <- RNA_4[-testIndexes, ] # Rest of the data 
  fit2 <- lm(fmla, data = trainData) # Trained using same formula but now on 4/5th of the data
  predpr<- predict(fit2,newdata =testData, type=c("response")) # Prediction in test data 
  testData$prob2=predpr
  t = cor.test(testData$prob2, testData$FOLD)
  print(ggplot(testData, aes(FOLD,prob2)) + 
          geom_point() + 
          ggtitle (paste(i,": R=", round(t$estimate, digits = 2))) + 
          xlab("Measured Relative RNA Expression (log2)") + 
          ylab("Predicted Relative RNA Expression (log2)") + 
          geom_smooth(method = lm) +   geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
          geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
          theme_science())
  ggsave(paste("cross_validation_all_interactions_RNA",i,".tiff" ,sep =""), width=5, height = 5, unit = "in")
}

RNA_4$residual <- RNA_4$FOLD - RNA_4$problm

#plot residuals of RNA model
ggplot(RNA_4, aes(problm,residual)) + 
  geom_point() + 
  xlab("Predicted Relative RNA Expression") + 
  ylab("Residual") + 
  geom_smooth(method = lm) +   geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()
ggsave("RNA_residual_plot.tiff", height=4, width=4, units = "in")


#---Re-format Inputs------------------------------------------------------------------------------------
# Calculate mean FOLD expression for all barcodes for a given pattern
aggregate(FOLD ~ RE_Identity, data= TE, FUN=median) -> TE_median
aggregate(FOLD ~ RE_Identity, data= TE, FUN=mean) -> TE_mean
# Add Sequence information to pattern  
TE_median = merge(TE_median, Seq, by ="RE_Identity", all.x = TRUE)
TE = TE_median
# Change control to BBBB
TE[,1] <- as.character(TE[,1])
TE[which(TE[,1] == "Control"),1] <- ("BBBB")

colnames(TE) <- c("RE_Identity","FOLD","sequence")

# Split the pattern to get individual positions 
X =(strsplit(as.character(TE$RE_Identity), split= ""))
n.obs <- sapply(X, length) # deals with different column numbers in each row
seq.max <- seq_len(max(n.obs))# deals with different column numbers in each row
X <- t(sapply(X, "[", i = seq.max))# deals with different column numbers in each row
X[X == 7] <- "L" # Replace 7 with L 
colnames(X) = paste("P", c(1:8), sep ="") #Change column names 
TE = cbind(TE, X) 
TE$type = nchar((as.character(TE$RE_Identity)))
TE_4 = subset(TE, type == 4) # Keep guys will 4 position pattern only
TE_4 = subset(TE_4, RE_Identity != "Dna2" ) # Remove Dna2 guy
#TE_4<- subset(TE_4, !grepl("h", TE_4$RE_Identity))

#---Linear Models-----------------------------------------------------------------------------------
# Set Reference level to B for linear model fits 
TE_4 <- within(TE_4, P1 <- relevel(P1, ref = "B"))
TE_4 <- within(TE_4, P2 <- relevel(P2, ref = "B"))
TE_4 <- within(TE_4, P3 <- relevel(P3, ref = "B"))
TE_4 <- within(TE_4, P4 <- relevel(P4, ref = "B"))

# Formula for the linear model 
fmla <- as.formula(paste("FOLD ~ P1*P2 + P1*P3 + P1*P4 + P2*P3 + P2*P4 + P3*P4  "))
# this particular model tries to predict expression with all 4 positions and all interactions
# P1*P2 + P1*P3 + P1*P4 + P2*P3 + P2*P4 + P3*P4
# Fit all data
fit <-  glm(fmla, data = TE_4)
summary(fit) #displays summary to the fit. Information stored in this variable
predpr<- predict(fit,type=c("response")) # Predicted values from the fit 
TE_4$problm=predpr # Adding predicted values to the dataset
t= cor(TE_4$problm, TE_4$FOLD)
ggplot(TE_4, aes(FOLD,problm)) + geom_point() + 
  ggtitle(label=(paste("R =", format(t, digits = 2)))) +
  xlab("Measured Relative TE (log2)") + ylab("Predicted Relative TE (log2)") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()
ggsave("Fit_all_interactions_TE.tiff", width=5, height = 5, unit = "in")

#calculate Cohen's F2
f2_TE= (t^2)/(1-t^2)
#perform power analysis
pwr.f2.test(u=112, v=512, f2= f2_TE, sig.level = 0.01, power=NULL)

data.frame(summary(fit)$coeff)-> coeff # stores information about paramters, estimates and significance 
coeff_TE=coeff[with(coeff, order(coeff[,4])), ] # you can order coeffs by Effect (Column1) or Signficance (Column 4)


# Perform fit 5 fold cross validation
TE_4<-TE_4[sample(nrow(TE_4)),]  # Randomize rows
folds <- cut(seq(1,nrow(TE_4)),breaks=5,labels=FALSE) #Assigns every row to be one part of 5 folds

for(i in 1:5){
  testIndexes <- which(folds==i,arr.ind=TRUE) #Segement your data by fold using the which() function 
  testData <-TE_4[testIndexes, ] # Data assinged to current index
  trainData <- TE_4[-testIndexes, ] # Rest of the data 
  fit2 <- lm(fmla, data = trainData) # Trained using same formula but now on 4/5th of the data
  predpr<- predict(fit2,newdata =testData, type=c("response")) # Prediction in test data 
  testData$prob2=predpr
  t = cor.test(testData$prob2, testData$FOLD)
  print(ggplot(testData, aes(FOLD,prob2)) + 
          geom_point() + 
          ggtitle (paste(i,": R=", round(t$estimate, digits = 2))) + 
          xlab("Measured Relative TE (log2)") + 
          ylab("Predicted Relative TE (log2)") + 
          geom_smooth(method = lm) +   geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
          geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
          theme_science())
  ggsave(paste("cross_validation_all_interactions_TE",i,".tiff" ,sep =""), width=5, height = 5, unit = "in")
}

#calculate residuals
TE_4$residual <- TE_4$FOLD - TE_4$problm
#residual plot
ggplot(TE_4, aes(problm,residual)) + 
  geom_point() + 
  xlab("Predicted Relative TE") + 
  ylab("Residual") + 
  geom_smooth(method = lm) +   geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()
ggsave("TE_residual_plot.tiff", height=4, width=4, units = "in")

#model predictions for onyl let-7 targets filter TE_4 and find only Let-7 targets
let_7_TE <- dplyr::filter(TE_4, RE_Identity %in% bcid_rename$Let.7)
#correlation between predicted v observed
t_let_te= cor.test(let_7_TE$problm, let_7_TE$FOLD)
#plot predicted v observed
ggplot(let_7_TE, aes(FOLD,problm)) + geom_point() + 
  ggtitle (paste("R=", round(t_let_te$estimate, digits = 3))) + 
  xlab("Measured Relative TE (log2)") + 
  ylab("Predicted Relative TE (log2)") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  theme_science()
ggsave("model_prediction_let_TE.tiff", height=6, width=6, units="in")

#model predictions for onyl let-7 targets filter TE_4 and find only Pum targets find correlation
pum_TE <- dplyr::filter(TE_4, RE_Identity %in% bcid_rename$Pumilio)
t_pum_te= cor.test(pum_TE$problm, pum_TE$FOLD)

#model predictions for onyl let-7 targets filter TE_4 and find only ARE targets find correlation
hur_TE <- dplyr::filter(TE_4, RE_Identity %in% bcid_rename$HuR)
t_hur_te= cor.test(hur_TE$problm, hur_TE$FOLD)

#model predictions for onyl let-7 targets filter TE_4 and find only Smg targets find correlation
smg_TE <- dplyr::filter(TE_4, RE_Identity %in% bcid_rename$Smg)
t_smg_te= cor.test(smg_TE$problm, smg_TE$FOLD)

#model predictions for onyl let-7 targets filter TE_4 and find only Pum-Let targets find correlation
pum_let_TE <- dplyr::filter(TE_4, RE_Identity %in% bcid_rename$Pum_Let_7)
t_pum_let_te= cor.test(pum_let_TE$problm, pum_let_TE$FOLD)

#model predictions for onyl let-7 targets filter TE_4 and find only Pum-ARE targets find correlation
pum_hur_TE <- dplyr::filter(TE_4, RE_Identity %in% bcid_rename$HuR_Pum)
t_pum_hur_te= cor.test(pum_hur_TE$problm, pum_hur_TE$FOLD)

#model predictions for onyl let-7 targets filter TE_4 and find only HuR-ARE targets find correlation
let_hur_TE <- dplyr::filter(TE_4, RE_Identity %in% bcid_rename$HuR_Let_7)
t_let_hur_te= cor.test(let_hur_TE$problm, let_hur_TE$FOLD)

#Make a summary table for all correlations for predicted v observed for each group
summary_stats <- data.frame(Regulatory_Element <- c("Let-7", "PRE", "ARE", "SRE", "Let-7/PRE", "Let-7/ARE", "ARE/PRE"), 
                            Correlation_RNA <- c(t_let$estimate, t_pum$estimate, t_hur$estimate, t_smg$estimate, t_pum_let$estimate, t_let_hur$estimate, t_pum_hur$estimate),
                            Correlation_TE <- c(t_let_te$estimate, t_pum_te$estimate, t_hur_te$estimate, t_smg_te$estimate, t_pum_let_te$estimate, t_let_hur_te$estimate, t_pum_hur_te$estimate))

write.csv(summary_stats, "summary_stats.csv")



#Make new column in coeff data.frames for ID
coeff_RNA$ID <- rownames(coeff_RNA)
coeff_TE$ID <- rownames(coeff_TE)


#merge and RNA and TE coeff 
coeff_RNA_TE <- merge(coeff_RNA, coeff_TE, by="ID")
colnames(coeff_RNA_TE) <- c("ID", "Estimate_RNA","Std.Error_RNA","t.value_RNA", "p.value_RNA", "Estimate_TE","Std.Error_TE","t.value_TE", "p.value_TE")

#correlation between RNA and TE coeff.
r=cor(coeff_RNA_TE$Estimate_RNA, coeff_RNA_TE$Estimate_TE, method="pearson")
#plot RNA v TE coeff.
ggplot(coeff_RNA_TE, aes(Estimate_RNA, Estimate_TE)) + 
  geom_point() + 
  # scale_x_continuous(limits=c(-1.5,1.8), breaks=c(-1, 0, 1)) +
  # scale_y_continuous(limits=c(-1.5,1.8), breaks=c(-1, 0, 1)) +
  geom_smooth(method="lm", colour="#56B4E9", se = FALSE) +
  geom_smooth(colour="#E69F00", se = FALSE) +
  xlab("RNA Coeff.") + 
  ylab("TE Coeff.") +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") +
  #geom_abline() +
  annotate("text", x=-1, y=0.3, fontface=2, label=(paste("R =", format(r, digits = 2)))) +
  theme_science()

ggsave("coeff_RNA_TE.tiff", width=5, height=4, units="in")

#subset coeff that are significant
coeff_RNA_sig <- subset(coeff_RNA, Pr...t.. < 0.05)
coeff_TE_sig <- subset(coeff_TE, Pr...t.. < 0.05)

#Assign level to coeff for RNA or TE
coeff_RNA$Level <- "RNA"
coeff_TE$Level <- "TE"
#bind coeff for RNA and TE
coeff_RNA_TE <- rbind(coeff_RNA, coeff_TE)
#remove intercept
coeff_RNA_TE <- subset(coeff_RNA_TE, !ID == "(Intercept)")
#rename h -> A
coeff_RNA_TE$ID <- gsub("h", "A", coeff_RNA_TE$ID)
#write table for all coeff from RNA and TE
write.csv(coeff_RNA_TE, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_TE.csv", row.names=FALSE)
#subset coeff for only let-7 targets
coeff_RNA_let <- subset(coeff_RNA_TE, !grepl("P.S|P.p|P.A", coeff_RNA_TE$ID))
write.csv(coeff_RNA_let, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_let.csv", row.names=FALSE)
#subset coeff for only pum targets
coeff_RNA_pum <- subset(coeff_RNA_TE, !grepl("P.S|P.L|P.A", coeff_RNA_TE$ID))
write.csv(coeff_RNA_pum, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_pum.csv", row.names=FALSE)
#subset coeff for only smg targets
coeff_RNA_smg <- subset(coeff_RNA_TE, !grepl("P.L|P.p|P.A", coeff_RNA_TE$ID))
write.csv(coeff_RNA_smg, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_smg.csv", row.names=FALSE)
#subset coeff for only ARE targets
coeff_RNA_hur <- subset(coeff_RNA_TE, !grepl("P.S|P.p|P.L", coeff_RNA_TE$ID))
write.csv(coeff_RNA_hur, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_hur.csv", row.names=FALSE)
#subset coeff for only Let-Pum targets
coeff_RNA_let_pum <- subset(coeff_RNA_TE, !grepl("P.S|P.A|P.L[:]P.L|P.p[:]P.p", coeff_RNA_TE$ID))
write.csv(coeff_RNA_let_pum, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_let_pum.csv", row.names=FALSE)
#subset coeff for only Let-Smg targets
coeff_RNA_let_smg <- subset(coeff_RNA_TE, !grepl("P.p|P.A|P.L[:]P.L|P.S[:]P.S", coeff_RNA_TE$ID))
write.csv(coeff_RNA_let_smg, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_let_smg.csv", row.names=FALSE)
#subset coeff for only Let-ARE targets
coeff_RNA_let_hur <- subset(coeff_RNA_TE, !grepl("P.S|P.p|P.L[:]P.L|P.A[:]P.A", coeff_RNA_TE$ID))
write.csv(coeff_RNA_let_hur, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_let_hur.csv", row.names=FALSE)
#subset coeff for only Pum-ARE targets
coeff_RNA_pum_hur <- subset(coeff_RNA_TE, !grepl("P.S|P.L|P.A[:]P.A|P.p[:]P.p", coeff_RNA_TE$ID))
write.csv(coeff_RNA_pum_hur, "/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_pum_hur.csv", row.names=FALSE)

#load packages
library(scales)
library(reshape)
library(plyr)

theme_science2 <- function (base_size = 18, base_family = "Arial Black") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=2), 
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(),
          axis.line.x = element_blank(),  axis.line.y = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}



#load cbpalette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#make a new column labeld "Sig" that contains p.values
coeff_RNA_let_pum$Sig <- coeff_RNA_let_pum$Pr...t.. 
#covert p.value in to *, **, or ***
coeff_RNA_let_pum$Sig[coeff_RNA_let_pum$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_let_pum$Sig[coeff_RNA_let_pum$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_let_pum$Sig[coeff_RNA_let_pum$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_let_pum$Sig[coeff_RNA_let_pum$Pr...t.. > 0.05] <- as.character("")

#plot pum-let coeff for RNA and TE
ggplot(coeff_RNA_let_pum, aes(reorder(ID, Estimate), Level)) + 
  geom_tile(aes(fill = Estimate), colour="white") + 
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.") +
  xlab("") + ylab("") +
  annotate("text", x=coeff_RNA_let_pum$ID[coeff_RNA_let_pum$Level == "RNA"], y=1, fontface=4, label=(paste("", coeff_RNA_let_pum$Sig[coeff_RNA_let_pum$Level == "RNA"]))) +
  annotate("text", x=coeff_RNA_let_pum$ID[coeff_RNA_let_pum$Level == "TE"], y=2, fontface=4, label=(paste("", coeff_RNA_let_pum$Sig[coeff_RNA_let_pum$Level == "TE"]))) +
  coord_flip() +
  theme_science2()
ggsave("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/let_pum_interactions.tiff", width=5, height=5, units="in")

#make a new column labeld "Sig" that contains p.values
coeff_RNA_let_hur$Sig <- coeff_RNA_let_hur$Pr...t.. 
#covert p.value in to *, **, or ***
coeff_RNA_let_hur$Sig[coeff_RNA_let_hur$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_let_hur$Sig[coeff_RNA_let_hur$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_let_hur$Sig[coeff_RNA_let_hur$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_let_hur$Sig[coeff_RNA_let_hur$Pr...t.. > 0.05] <- as.character("")
#plot ARE-let coeff for RNA and TE
ggplot(coeff_RNA_let_hur, aes(reorder(ID, Estimate), Level)) + 
  geom_tile(aes(fill = Estimate), colour="white") + 
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.") +
  xlab("") + ylab("") +
  annotate("text", x=coeff_RNA_let_hur$ID[coeff_RNA_let_hur$Level == "RNA"], y=1, fontface=4, label=(paste("", coeff_RNA_let_hur$Sig[coeff_RNA_let_hur$Level == "RNA"]))) +
  annotate("text", x=coeff_RNA_let_hur$ID[coeff_RNA_let_hur$Level == "TE"], y=2, fontface=4, label=(paste("", coeff_RNA_let_hur$Sig[coeff_RNA_let_hur$Level == "TE"]))) +
  coord_flip() +
  theme_science2()
ggsave("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/let_hur_interactions.tiff", width=5, height=5, units="in")

#make a new column labeld "Sig" that contains p.values
coeff_RNA_pum_hur$Sig <- coeff_RNA_pum_hur$Pr...t.. 
#covert p.value in to *, **, or ***
coeff_RNA_pum_hur$Sig[coeff_RNA_pum_hur$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_pum_hur$Sig[coeff_RNA_pum_hur$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_pum_hur$Sig[coeff_RNA_pum_hur$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_pum_hur$Sig[coeff_RNA_pum_hur$Pr...t.. > 0.05] <- as.character("")
#plot ARE-pum coeff for RNA and TE
ggplot(coeff_RNA_pum_hur, aes(reorder(ID, Estimate), Level)) + 
  geom_tile(aes(fill = Estimate), colour="white") + 
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.") +
  xlab("") + ylab("") +
  annotate("text", x=coeff_RNA_pum_hur$ID[coeff_RNA_pum_hur$Level == "RNA"], y=1, fontface=4, label=(paste("", coeff_RNA_pum_hur$Sig[coeff_RNA_pum_hur$Level == "RNA"]))) +
  annotate("text", x=coeff_RNA_pum_hur$ID[coeff_RNA_pum_hur$Level == "TE"], y=2, fontface=4, label=(paste("", coeff_RNA_pum_hur$Sig[coeff_RNA_pum_hur$Level == "TE"]))) +
  coord_flip() +
  theme_science2()
ggsave("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/pum_hur_interactions.tiff", width=5, height=5, units="in")

#make a new column labeld "Sig" that contains p.values
coeff_RNA_pum$Sig <- coeff_RNA_pum$Pr...t.. 
#covert p.value in to *, **, or ***
coeff_RNA_pum$Sig[coeff_RNA_pum$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_pum$Sig[coeff_RNA_pum$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_pum$Sig[coeff_RNA_pum$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_pum$Sig[coeff_RNA_pum$Pr...t.. > 0.05] <- as.character("")
#plot pum coeff for RNA and TE
ggplot(coeff_RNA_pum, aes(reorder(ID, Estimate), Level)) + 
  geom_tile(aes(fill = Estimate), colour="white") + 
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.") +
  xlab("") + ylab("") +
  annotate("text", x=coeff_RNA_pum$ID[coeff_RNA_pum$Level == "RNA"], y=1, fontface=4, label=(paste("", coeff_RNA_pum$Sig[coeff_RNA_pum$Level == "RNA"]))) +
  annotate("text", x=coeff_RNA_pum$ID[coeff_RNA_pum$Level == "TE"], y=2, fontface=4, label=(paste("", coeff_RNA_pum$Sig[coeff_RNA_pum$Level == "TE"]))) +
  coord_flip() +
  theme_science2()
ggsave("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/pum_interactions.tiff", width=5, height=5, units="in")

#make a new column labeld "Sig" that contains p.values
coeff_RNA_let$Sig <- coeff_RNA_let$Pr...t.. 
#covert p.value in to *, **, or ***
coeff_RNA_let$Sig[coeff_RNA_let$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_let$Sig[coeff_RNA_let$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_let$Sig[coeff_RNA_let$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_let$Sig[coeff_RNA_let$Pr...t.. > 0.05] <- as.character("")
#plot let coeff for RNA and TE
ggplot(coeff_RNA_let, aes(reorder(ID, Estimate), Level)) + 
  geom_tile(aes(fill = Estimate), colour="white") + 
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.") +
  xlab("") + ylab("") +
  annotate("text", x=coeff_RNA_let$ID[coeff_RNA_let$Level == "RNA"], y=1, fontface=4, label=(paste("", coeff_RNA_let$Sig[coeff_RNA_let$Level == "RNA"]))) +
  annotate("text", x=coeff_RNA_let$ID[coeff_RNA_let$Level == "TE"], y=2, fontface=4, label=(paste("", coeff_RNA_let$Sig[coeff_RNA_let$Level == "TE"]))) +
  coord_flip() +
  theme_science2()
ggsave("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/let_interactions.tiff", width=5, height=5, units="in")

#make a new column labeld "Sig" that contains p.values
coeff_RNA_hur$Sig <- coeff_RNA_hur$Pr...t.. 
#covert p.value in to *, **, or ***
coeff_RNA_hur$Sig[coeff_RNA_hur$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_hur$Sig[coeff_RNA_hur$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_hur$Sig[coeff_RNA_hur$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_hur$Sig[coeff_RNA_hur$Pr...t.. > 0.05] <- as.character("")
#plot ARE coeff for RNA and TE
ggplot(coeff_RNA_hur, aes(reorder(ID, Estimate), Level)) + 
  geom_tile(aes(fill = Estimate), colour="white") + 
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.") +
  xlab("") + ylab("") +
  annotate("text", x=coeff_RNA_hur$ID[coeff_RNA_hur$Level == "RNA"], y=1, fontface=4, label=(paste("", coeff_RNA_hur$Sig[coeff_RNA_hur$Level == "RNA"]))) +
  annotate("text", x=coeff_RNA_hur$ID[coeff_RNA_hur$Level == "TE"], y=2, fontface=4, label=(paste("", coeff_RNA_hur$Sig[coeff_RNA_hur$Level == "TE"]))) +
  coord_flip() +
  theme_science2()
ggsave("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/hur_interactions.tiff", width=5, height=5, units="in")

#read coeff for let-pum annotated (annotated by hand)
coeff_RNA_let_pum_comb <- read.csv("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_let_pum_annotated.csv")
#covert p.value in to *, **, or ***
coeff_RNA_let_pum_comb$Sig <- coeff_RNA_let_pum_comb$Pr...t.. 
coeff_RNA_let_pum_comb$Sig[coeff_RNA_let_pum_comb$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_let_pum_comb$Sig[coeff_RNA_let_pum_comb$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_let_pum_comb$Sig[coeff_RNA_let_pum_comb$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_let_pum_comb$Sig[coeff_RNA_let_pum_comb$Pr...t.. > 0.05] <- as.character("")


#plot RNA coeff for let-pum as matrix
ggplot(subset(coeff_RNA_let_pum_comb, Level == "RNA"), aes(PRE, LET, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.\nRNA") +
  xlab("PRE Position") + ylab("Let-7 Position")  +
  annotate("text", x=coeff_RNA_let_pum_comb$PRE[coeff_RNA_let_pum_comb$Level == "RNA"], y=coeff_RNA_let_pum_comb$LET[coeff_RNA_let_pum_comb$Level == "RNA"], fontface=4, label=(paste("", coeff_RNA_let_pum_comb$Sig[coeff_RNA_let_pum_comb$Level == "RNA"]))) +
  theme_science2()
ggsave("let_pum_coeff_RNA_matrix.tiff", height=4, width =5, units="in")

#plot TE coeff for let-pum as matrix
ggplot(subset(coeff_RNA_let_pum_comb, Level == "TE"), aes(PRE, LET, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.\nTE") +
  xlab("PRE Position") + ylab("Let-7 Position") +
  annotate("text", x=coeff_RNA_let_pum_comb$PRE[coeff_RNA_let_pum_comb$Level == "TE"], y=coeff_RNA_let_pum_comb$LET[coeff_RNA_let_pum_comb$Level == "TE"], fontface=4, label=(paste("", coeff_RNA_let_pum_comb$Sig[coeff_RNA_let_pum_comb$Level == "TE"]))) +
  theme_science2()
ggsave("let_pum_coeff_TE_matrix.tiff", height=4, width =5, units="in")

#read coeff for ARE (HuR)-let annotated (annotated by hand)
coeff_RNA_let_hur_comb <- read.csv("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_let_hur_annotated.csv")
#covert p.value in to *, **, or ***
coeff_RNA_let_hur_comb$Sig <- coeff_RNA_let_hur_comb$Pr...t.. 
coeff_RNA_let_hur_comb$Sig[coeff_RNA_let_hur_comb$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_let_hur_comb$Sig[coeff_RNA_let_hur_comb$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_let_hur_comb$Sig[coeff_RNA_let_hur_comb$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_let_hur_comb$Sig[coeff_RNA_let_hur_comb$Pr...t.. > 0.05] <- as.character("")


#plot RNA coeff for ARE-let as matrix
ggplot(subset(coeff_RNA_let_hur_comb, Level == "RNA"), aes(ARE, LET, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.\nRNA") +
  xlab("ARE Position") + ylab("Let-7 Position")  +
  annotate("text", x=coeff_RNA_let_hur_comb$ARE[coeff_RNA_let_hur_comb$Level == "RNA"], y=coeff_RNA_let_hur_comb$LET[coeff_RNA_let_hur_comb$Level == "RNA"], fontface=4, label=(paste("", coeff_RNA_let_hur_comb$Sig[coeff_RNA_let_hur_comb$Level == "RNA"]))) +
  theme_science2()
ggsave("let_hur_coeff_RNA_matrix.tiff", height=4, width =5, units="in")

#plot TE coeff for ARE-let as matrix
ggplot(subset(coeff_RNA_let_hur_comb, Level == "TE"), aes(ARE, LET, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.\nTE") +
  xlab("ARE Position") + ylab("Let-7 Position") +
  annotate("text", x=coeff_RNA_let_hur_comb$ARE[coeff_RNA_let_hur_comb$Level == "TE"], y=coeff_RNA_let_hur_comb$LET[coeff_RNA_let_hur_comb$Level == "TE"], fontface=4, label=(paste("", coeff_RNA_let_hur_comb$Sig[coeff_RNA_let_hur_comb$Level == "TE"]))) +
  theme_science2()       
ggsave("let_hur_coeff_TE_matrix.tiff", height=4, width =5, units="in")

#read coeff for ARE (HuR)-pum annotated (annotated by hand)
coeff_RNA_pum_hur_comb <- read.csv("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_pum_hur_annotated.csv")
#covert p.value in to *, **, or ***
coeff_RNA_pum_hur_comb$Sig <- coeff_RNA_pum_hur_comb$Pr...t.. 
coeff_RNA_pum_hur_comb$Sig[coeff_RNA_pum_hur_comb$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_pum_hur_comb$Sig[coeff_RNA_pum_hur_comb$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_pum_hur_comb$Sig[coeff_RNA_pum_hur_comb$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_pum_hur_comb$Sig[coeff_RNA_pum_hur_comb$Pr...t.. > 0.05] <- as.character("")

#plot RNA coeff for ARE-pum as matrix
ggplot(subset(coeff_RNA_pum_hur_comb, Level == "RNA"), aes(ARE, PRE, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.\nRNA") +
  xlab("ARE Position") + ylab("PRE Position")  +
  annotate("text", x=coeff_RNA_pum_hur_comb$ARE[coeff_RNA_pum_hur_comb$Level == "RNA"], y=coeff_RNA_pum_hur_comb$PRE[coeff_RNA_pum_hur_comb$Level == "RNA"], fontface=4, label=(paste("", coeff_RNA_pum_hur_comb$Sig[coeff_RNA_pum_hur_comb$Level == "RNA"]))) +
  theme_science2()
ggsave("pum_hur_coeff_RNA_matrix.tiff", height=4, width =5, units="in")

#plot TE coeff for ARE-pum as matrix
ggplot(subset(coeff_RNA_pum_hur_comb, Level == "TE"), aes(ARE, PRE, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff.\nTE") +
  xlab("ARE Position") + ylab("PRE Position") +
  annotate("text", x=coeff_RNA_pum_hur_comb$ARE[coeff_RNA_pum_hur_comb$Level == "TE"], y=coeff_RNA_pum_hur_comb$PRE[coeff_RNA_pum_hur_comb$Level == "TE"], fontface=4, label=(paste("", coeff_RNA_pum_hur_comb$Sig[coeff_RNA_pum_hur_comb$Level == "TE"]))) +
  theme_science2()    
ggsave("pum_hur_coeff_TE_matrix.tiff", height=4, width =5, units="in")

#read coeff for let-7 annotated (annotated by hand)
coeff_RNA_let_comb <- read.csv("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_let_annotated.csv")
#covert p.value in to *, **, or ***
coeff_RNA_let_comb$Sig <- coeff_RNA_let_comb$Pr...t.. 
coeff_RNA_let_comb$Sig[coeff_RNA_let_comb$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_let_comb$Sig[coeff_RNA_let_comb$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_let_comb$Sig[coeff_RNA_let_comb$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_let_comb$Sig[coeff_RNA_let_comb$Pr...t.. > 0.05] <- as.character("")

#plot RNA coeff for let-7 as matrix
ggplot(subset(coeff_RNA_let_comb, Level == "RNA"), aes(Let.2, Let.1, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff. RNA") +
  xlab("Let-7 Position") + ylab("Let-7 Position")  +
  annotate("text", x=coeff_RNA_let_comb$Let.2[coeff_RNA_let_comb$Level == "RNA"], y=coeff_RNA_let_comb$Let.1[coeff_RNA_let_comb$Level == "RNA"], fontface=4, angle=45, label=(paste("", coeff_RNA_let_comb$Sig[coeff_RNA_let_comb$Level == "RNA"]))) +
  theme_science2() + theme(legend.justification = c(1, 0), legend.position = c(0.6, 0.75), legend.direction = "horizontal") + guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5))
ggsave("let_coeff_RNA_matrix.tiff", height=4, width =4, units="in")

#plot TE coeff for let-7 as matrix
ggplot(subset(coeff_RNA_let_comb, Level == "RNA"), aes(Let.2, Let.1, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff. TE") +
  xlab("Let-7 Position") + ylab("Let-7 Position")  +
  annotate("text", x=coeff_RNA_let_comb$Let.2[coeff_RNA_let_comb$Level == "TE"], y=coeff_RNA_let_comb$Let.1[coeff_RNA_let_comb$Level == "TE"], fontface=4, angle=45, label=(paste("", coeff_RNA_let_comb$Sig[coeff_RNA_let_comb$Level == "TE"]))) +
  theme_science2() + theme(legend.justification = c(1, 0), legend.position = c(0.6, 0.75), legend.direction = "horizontal") + guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5))
ggsave("let_coeff_TE_matrix.tiff", height=4, width =4, units="in")

#read coeff for pum annotated (annotated by hand)
coeff_RNA_pum_comb <- read.csv("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_pum_annotated.csv")
#covert p.value in to *, **, or ***
coeff_RNA_pum_comb$Sig <- coeff_RNA_pum_comb$Pr...t.. 
coeff_RNA_pum_comb$Sig[coeff_RNA_pum_comb$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_pum_comb$Sig[coeff_RNA_pum_comb$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_pum_comb$Sig[coeff_RNA_pum_comb$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_pum_comb$Sig[coeff_RNA_pum_comb$Pr...t.. > 0.05] <- as.character("")

#plot RNA coeff for pum as matrix
ggplot(subset(coeff_RNA_pum_comb, Level == "RNA"), aes(Pum.2, Pum.1, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff. RNA") +
  xlab("PRE Position") + ylab("PRE Position")  +
  annotate("text", x=coeff_RNA_pum_comb$Pum.2[coeff_RNA_pum_comb$Level == "RNA"], y=coeff_RNA_pum_comb$Pum.1[coeff_RNA_pum_comb$Level == "RNA"], fontface=4, angle=45, label=(paste("", coeff_RNA_pum_comb$Sig[coeff_RNA_pum_comb$Level == "RNA"]))) +
  theme_science2() + theme(legend.justification = c(1, 0), legend.position = c(0.6, 0.75), legend.direction = "horizontal") + guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5))
ggsave("pum_coeff_RNA_matrix.tiff", height=4, width =4, units="in")

#plot te coeff for pum as matrix
ggplot(subset(coeff_RNA_pum_comb, Level == "RNA"), aes(Pum.2, Pum.1, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff. TE") +
  xlab("PRE Position") + ylab("PRE Position")  +
  annotate("text", x=coeff_RNA_pum_comb$Pum.2[coeff_RNA_pum_comb$Level == "TE"], y=coeff_RNA_pum_comb$Pum.1[coeff_RNA_pum_comb$Level == "TE"], fontface=4, angle=45, label=(paste("", coeff_RNA_pum_comb$Sig[coeff_RNA_pum_comb$Level == "TE"]))) +
  theme_science2() + theme(legend.justification = c(1, 0), legend.position = c(0.6, 0.75), legend.direction = "horizontal") + guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5))
ggsave("pum_coeff_TE_matrix.tiff", height=4, width =4, units="in")

#read coeff for ARE (HuR) annotated (annotated by hand)
coeff_RNA_hur_comb <- read.csv("/Users/Kyle/Dropbox/CRE_seq/RE_Array_EX1/coeff_RNA_hur_annotated.csv")
#covert p.value in to *, **, or ***
coeff_RNA_hur_comb$Sig <- coeff_RNA_hur_comb$Pr...t.. 
coeff_RNA_hur_comb$Sig[coeff_RNA_hur_comb$Pr...t.. < 0.05] <- as.character("*")
coeff_RNA_hur_comb$Sig[coeff_RNA_hur_comb$Pr...t.. < 0.01] <- as.character("**")
coeff_RNA_hur_comb$Sig[coeff_RNA_hur_comb$Pr...t.. < 0.001] <- as.character("***")
coeff_RNA_hur_comb$Sig[coeff_RNA_hur_comb$Pr...t.. > 0.05] <- as.character("")

#plot RNA coeff for ARE as matrix
ggplot(subset(coeff_RNA_hur_comb, Level == "RNA"), aes(hur.2, hur.1, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff. RNA") +
  xlab("ARE Position") + ylab("ARE Position")  +
  annotate("text", x=coeff_RNA_hur_comb$hur.2[coeff_RNA_hur_comb$Level == "RNA"], y=coeff_RNA_hur_comb$hur.1[coeff_RNA_hur_comb$Level == "RNA"], fontface=4, angle=45, label=(paste("", coeff_RNA_hur_comb$Sig[coeff_RNA_hur_comb$Level == "RNA"]))) +
  theme_science2() + theme(legend.justification = c(1, 0), legend.position = c(0.6, 0.75), legend.direction = "horizontal") + guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5))
ggsave("hur_coeff_RNA_matrix.tiff", height=4, width =4, units="in")

#plot te coeff for ARE as matrix
ggplot(subset(coeff_RNA_hur_comb, Level == "RNA"), aes(hur.2, hur.1, fill=Estimate)) + geom_tile()  +
  scale_fill_gradient2(low ="#E69F00", high =  "#56B4E9", name = "Coeff. TE") +
  xlab("ARE Position") + ylab("ARE Position")  +
  annotate("text", x=coeff_RNA_hur_comb$hur.2[coeff_RNA_hur_comb$Level == "TE"], y=coeff_RNA_hur_comb$hur.1[coeff_RNA_hur_comb$Level == "TE"], fontface=4, angle=45, label=(paste("", coeff_RNA_hur_comb$Sig[coeff_RNA_hur_comb$Level == "TE"]))) +
  theme_science2() + theme(legend.justification = c(1, 0), legend.position = c(0.6, 0.75), legend.direction = "horizontal") + guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5))
ggsave("hur_coeff_TE_matrix.tiff", height=4, width =4, units="in")

#load thermo model for RNA outputs
thermo_RNA <- read.delim("Kyle_thermo_output_exp")
#plot expected v observed for thermo model of RNA
ggplot(thermo_RNA, aes(Observed, Predicted)) + geom_point() + theme_science() + 
  xlab("Measured Relative RNA Expression (log2)") + 
  ylab("Predicted Relative RNA Expression (log2)") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") 
ggsave("thermo_RNA.tiff", height=4.5, width=4.5, units="in", dpi=300)

#load thermo model for TE outputs
thermo_TE <- read.delim("Kyle_thermo_output_exp_TE")
#plot expected v observed for thermo model of TE
ggplot(thermo_RNA, aes(Observed, Predicted)) + geom_point() + theme_science() + 
  xlab("Measured Relative TE (log2)") + 
  ylab("Predicted Relative TE (log2)") + 
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, linetype="dashed", colour = "gray") +
  geom_vline(xintercept = 0, linetype="dashed", colour = "gray") 
ggsave("thermo_TE.tiff", height=4.5, width=4.5, units="in", dpi=300)