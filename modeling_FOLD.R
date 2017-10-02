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
library("stringr")

#---Set Working Directory (Input Output Files)------------------------------------------------------------
setwd("/Users/Hemangi/Box Sync/2017/Kyle/")

#---INPUT-------------------------------------------------------------------------------------------------
bcid = read.delim("/Users/Hemangi/Box Sync/2017/Kyle/Grouped_Barcode_identities.txt")
TE = read.table("/Users/Hemangi/Box Sync/2017/Kyle/PARNA_RNA.txt",header = TRUE)
RNA = read.table("/Users/Hemangi/Box Sync/2017/Kyle/RNA_plasmid.txt",header = TRUE)
Seq = read.table("/Users/Hemangi/Box Sync/2017/Kyle/Library.txt",header = TRUE)
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
fmla <- as.formula(paste("FOLD ~ P1 + P2 + P3 + P4 + P1:P2 + P1:P3 + P1:P4 + P2:P3 + P2:P4 + P3:P4"))
fmla_all <- as.formula(paste("FOLD ~ P1 + P2 + P3 + P4 + P1:P2 + P1:P3 + P1:P4 + P2:P3 + P2:P4 + P3:P4 + P1:P2:P3 + P1:P2:P4 + P2:P3:P4 + P1:P3:P4"))
fmla_al <- as.formula(paste("FOLD ~ P1*P2*P3 + P1*P2*P4 + P4*P2*P3 "))

# Explore Model space #Feb 27,2017
#---All possible terms------------------------
x <- c("P1", "P2", "P3", "P4","P1*P2", "P1*P3", "P1*P4",  "P2*P3", "P2*P4", "P3*P4" )

#---Exploring all possible models-----------------------
# goal is to create all possible combinations of paramters. We have 10 terms total. 
# We can start with 1 term only - there are 10 options..
# for 10 terms - There is only one combination. We want to explore all combinations
# We also have to keep in mind that interaction parameters include singles by default. 
Model_space = NULL
for( i in seq(1:10)){ 
  combn(x,i) -> A  # creates all possible combinations of length i 
  for (j in seq(1:ncol(A))){ # for each combination j  of length i 
    fmla <- as.formula(paste("FOLD ~ ", paste(A[,j], collapse  = "+"))) # we create a new formula with j
    fit <-  lm(fmla, data = RNA_4) # we fit a new model with j
    Model_space = rbind(Model_space,  c(paste(A[,j], collapse  = "+"), summary(fit)$r.square))  # append model terms and R-squared value to model space 
  }
}
Model_space = data.frame(Model_space) 
Model_space= Model_space[with(Model_space, order(Model_space[,2])), ]
Model_space$numterms = str_count((Model_space$X1), pattern = "\\+")  # I am simply counting the number of plus-es for getting total number of terms. This is wrong becuase interactions add more terms. this is fixed later 
Model_space$numInt = str_count((Model_space$X1), pattern = "\\*") #  I am simply counting the number of stars for getting number of interacting terms.

X =(strsplit(as.character(Model_space$X1), split= "\\+"))
n.obs <- sapply(X, length) # deals with different column numbers in each row
seq.max <- seq_len(max(n.obs))# deals with different column numbers in each row
X <- t(sapply(X, "[", i = seq.max))# deals with different column numbers in each row
X = melt(X, id = rownames(X))
X1 = X[,c(1,3)]
X1$value <- as.factor(X1$value) # resets factors 
Y = dcast(X1, Var1~value)
Y[which(Y[,3] == 1), c(2,6)] <- 1  # If an interaction term is present, change the matrix such that individual values are included 
Y[which(Y[,4] == 1), c(2,9)] <- 1 # do it for all interactions 
Y[which(Y[,5] == 1), c(2,11)] <- 1
Y[which(Y[,7] == 1), c(6,9)] <- 1
Y[which(Y[,8] == 1), c(11,6)] <- 1
Y[which(Y[,10] == 1), c(11,9)] <- 1
Model_space = cbind(Model_space, Y)
Model_space$numterms = rowSums((Model_space[,c(6:15)]))
Model_space = Model_space[,c(2:14)]
Model_space = unique(Model_space) # keep only unique combinations 
Model_space1 = melt(Model_space[,c(1:2, 4:13)], id.vars=c("X2", "numterms"))
Model_space1 = subset(Model_space1, value == 1)

# Plots Number of Parameters in the model vs R-squared. Each dot it colored by the number of interacting terms. Each dot is one model 
ggplot(Model_space, aes(numterms, as.numeric(as.character(X2)), color= as.factor(numInt))) + geom_point(alpha = 0.5) + xlab("Total number of parameters in the model") + ylab("R-squared for fit model") + scale_color_discrete( guide = guide_legend(title = "Number of Interactions" ))

# Each panel shows number of interactions in the model. Within a panel, boxplots for each variable show the distributions of r-square values for models that include it. 
ggplot(Model_space1, aes(variable, as.numeric(as.character(X2)))) + geom_boxplot() + facet_wrap(~numterms) + ylab("R-squared for fit model") 

Model_space1= Model_space1[with(Model_space1, order(Model_space1[,1])), ]
# Each row is one of the 141 unique models. Each column in colored pink if that parameter was present in the model. 
ggplot(Model_space1, aes(variable, (as.character(X2)))) +  geom_tile(aes(fill = as.factor(value))) +  ylab("R-squared for fit model") + scale_fill_discrete( guide = guide_legend(title = "Present" )) + xlab("Parameter")
ggsave("Heatmap.pdf", width = 6, height = 12)

# this particular model tries to predict expression with all 4 positions and all interactions
# P1*P2 + P1*P3 + P1*P4 + P2*P3 + P2*P4 + P3*P4
# Fit all data
fit_all <-  lm(fmla_all, data = RNA_4)
fit_al <-  lm(fmla_al, data = RNA_4)
fit <-  lm(fmla, data = RNA_4)
AIC(fit, fit_al, fit_all)
library(MASS)
step <- stepAIC(fit_all, direction="both")
step$anova
summary(fit_all) #displays summary to the fit. Information stored in this variable
predpr<- predict(fit,type=c("response")) # Predicted values from the fit 
RNA_4$problm=predpr # Adding predicted values to the dataset
t= cor.test(log2(RNA_4$problm), log2(RNA_4$FOLD))
ggplot(RNA_4, aes(FOLD,problm)) + geom_point() + ggtitle (paste(i,": R=", round(t$estimate, digits = 2))) + xlab("Measured Relative RNA Expression") + ylab("Predicted Relative RNA Expression") + geom_smooth(method = lm)
ggsave("Fit_all_interactions_RNA.pdf")

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
  print(ggplot(testData, aes(FOLD,prob2)) + geom_point() + ggtitle (paste(i,": R=", round(t$estimate, digits = 2))) + xlab("Measured Relative RNA Expression") + ylab("Predicted Relative RNA Expression") + geom_smooth(method = lm))
  #ggsave(paste("cross_validation_all_interactions_RNA",i,".pdf" ,sep =""))
}

#---Remove HuR-----------------------------------------------------------------------------------
RNA_4<- subset(RNA_4, !grepl("h", RNA_4$RE_Identity))

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
t= cor.test(log2(RNA_4$problm), log2(RNA_4$FOLD))
ggplot(RNA_4, aes(FOLD,problm)) + geom_point() + ggtitle (paste(i,": R=", round(t$estimate, digits = 2))) + xlab("Measured Relative RNA Expression") + ylab("Predicted Relative RNA Expression") + geom_smooth(method = lm)
ggsave("Fit_all_interactions_RNA_noh.pdf")

data.frame(summary(fit)$coeff)-> coeff # stores information about paramters, estimates and significance 
coeff_RNA_noh=coeff[with(coeff, order(coeff[,4])), ] # you can order coeffs by Effect (Column1) or Signficance (Column 4)


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
  print(ggplot(testData, aes(FOLD,prob2)) + geom_point() + ggtitle (paste(i,": R=", round(t$estimate, digits = 2))) + xlab("Measured Relative RNA Expression") + ylab("Predicted Relative RNA Expression") + geom_smooth(method = lm))
  ggsave(paste("cross_validation_all_interactions_RNA_noh",i,".pdf" ,sep =""))
}

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
t= cor.test(log2(TE_4$problm), log2(TE_4$FOLD))
ggplot(TE_4, aes(FOLD,problm)) + geom_point() + ggtitle (paste(i,": R=", round(t$estimate, digits = 2))) + xlab("Measured Relative TE") + ylab("Predicted Relative TE") + geom_smooth(method = lm)
ggsave("Fit_all_interactions_TE.pdf")

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
  print(ggplot(testData, aes(FOLD,prob2)) + geom_point() + ggtitle (paste(i,": R=", round(t$estimate, digits = 2))) + xlab("Measured Relative TE") + ylab("Predicted Relative TE") + geom_smooth(method = lm))
  ggsave(paste("cross_validation_all_interactions_TE",i,".pdf" ,sep =""))
}


#---Remove HuR-----------------------------------------------------------------------------------
TE_4<- subset(TE_4, !grepl("h", TE_4$RE_Identity))

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
t= cor.test(log2(TE_4$problm), log2(TE_4$FOLD))
ggplot(TE_4, aes(FOLD,problm)) + geom_point() + ggtitle (paste(i,": R=", round(t$estimate, digits = 2))) + xlab("Measured Relative TE") + ylab("Predicted Relative TE") + geom_smooth(method = lm)
ggsave("Fit_all_interactions_TE_noh.pdf")

data.frame(summary(fit)$coeff)-> coeff # stores information about paramters, estimates and significance 
coeff_TE_noh=coeff[with(coeff, order(coeff[,4])), ] # you can order coeffs by Effect (Column1) or Signficance (Column 4)


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
  print(ggplot(testData, aes(FOLD,prob2)) + geom_point() + ggtitle (paste(i,": R=", round(t$estimate, digits = 2))) + xlab("Measured Relative TE") + ylab("Predicted Relative TE") + geom_smooth(method = lm))
  ggsave(paste("cross_validation_all_interactions_TE_noh",i,".pdf" ,sep =""))
}


coeff_RNA$ID <- rownames(coeff_RNA)
coeff_TE$ID <- rownames(coeff_TE)
coeff_RNA_noh$ID <- rownames(coeff_RNA_noh)
coeff_TE_noh$ID <- rownames(coeff_TE_noh)


coeff_RNA_TE <- merge(coeff_RNA, coeff_TE, by="ID")

ggplot(coeff_RNA_TE, aes(Estimate.x, Estimate.y)) + geom_point() + xlab("RNA Effect") + ylab("TE Effect")

cor.test(coeff_RNA_TE$Estimate.x, coeff_RNA_TE$Estimate.y)

coeff_RNA_sig <- subset(coeff_RNA, Pr...t.. < 0.05)
coeff_RNA_noh_sig <- subset(coeff_RNA_noh, Pr...t.. < 0.05)
coeff_TE_sig <- subset(coeff_TE, Pr...t.. < 0.05)
coeff_TE__noh_sig <- subset(coeff_TE_noh, Pr...t.. < 0.05)

coeff_RNA_let <- subset(coeff_RNA_sig, !grepl("P*S|P*p|P*h", coeff_RNA_sig$ID))

coeff_RNA_pum <- subset(coeff_RNA_sig, !grepl("P*S|P*L|P*h", coeff_RNA_sig$ID))

coeff_RNA_smg <- subset(coeff_RNA_sig, !grepl("P*L|P*p|P*h", coeff_RNA_sig$ID))

coeff_RNA_hur <- subset(coeff_RNA_sig, !grepl("P*S|P*p|P*L", coeff_RNA_sig$ID))

coeff_RNA_let_pum <- subset(coeff_RNA_sig, !grepl("P*S|P*h", coeff_RNA_sig$ID))

coeff_RNA_let_smg <- subset(coeff_RNA_sig, !grepl("P*p|P*h", coeff_RNA_sig$ID))

coeff_RNA_let_hur <- subset(coeff_RNA_sig, !grepl("P*p|P*S", coeff_RNA_sig$ID))

coeff_RNA_pum_hur <- subset(coeff_RNA_sig, !grepl("P*L|P*S", coeff_RNA_sig$ID))
