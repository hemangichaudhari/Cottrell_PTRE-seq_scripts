library(dplyr)
library(ggplot2)
HuR_Pum = read.delim("/Users/Kyle/Desktop/RE-Array_SIS2/HuR_Pum_Results.txt")

ggplot(HuR_Pum, aes(x = reorder(RE_Identity, -Fold), y= Fold)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + 
  scale_y_continuous(limits = c(0,8), breaks = 1:8, expand = c(0, 0)) +
  geom_bar(colour="black", fill="black", position=position_dodge(), stat="identity") + theme_classic() +
  theme(axis.title.x = element_text(colour = "black", size = 20, face = "bold"), 
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 45, hjust = 1,), 
        axis.text.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1))

HuR_Pum_Condensed <- dplyr::group_by(HuR_Pum, condensed)
HuR_Pum_Condensed <- dplyr::summarise(HuR_Pum_Condensed, FOLD=mean(Fold))

ggplot(HuR_Pum_Condensed, aes(x = reorder(condensed, -FOLD), y= FOLD)) + 
  xlab("RE Identity") + ylab("Relative Expression") +
  scale_x_discrete(expand = c(0.01, 0)) + 
  scale_y_continuous(limits = c(0,4), breaks = 1:4, expand = c(0, 0)) +
  geom_bar(colour="black", fill="black", position=position_dodge(), stat="identity") + theme_classic() +
  theme(axis.title.x = element_text(colour = "black", size = 20, face = "bold"), 
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(colour = "black", size = 15, face = "bold", angle = 45, hjust = 1,), 
        axis.text.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1))