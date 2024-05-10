#WRKY CUI analysis

#importing necessary libraries
library(tidyverse)
library(metan)
library(corrplot)
library(readxl)
library(writexl)
library(dplyr)

#reading codon usage dataset from codonw
cui_dataset <- read_excel("CUI.xlsx")
View(cui_dataset)

#reading codon usage dataset from caical
caical <- read_excel("CAICal_GC123.xlsx")
View(caical)

#merging and cleaning two datasets to get final dataset form further analysis

joined_df <- merge(cui_dataset,caical[c(1,25,39,53)],
                   by.x = 'title',by.y = 'sequences',
                   all.x = TRUE,all.y = TRUE)
view(joined_df)

#creating new columns with unique values 
joined_df$GC1 <- joined_df$`%G1+C1`/100
joined_df$GC2 <- joined_df$`%G2+C2`/100
joined_df$GC3 <- joined_df$`%G3+C3`/100
view(joined_df)

#generating final data frame for analysis
cui <- joined_df[-c(7:8,12:18)]
View(cui)
write.csv(cui,file = 'Supplimentary_datafile.csv',row.names = FALSE)

#checking normality of the variables by Shapiro-Wilk normality test

shapiro_wilk <- function(x){
  shapiro.test(x)$p.value} #creating function to perform Shapiro-Wilk test and return the p-values only
#
shapiro_result <- apply(cui[-c(1)],2,shapiro_wilk) #applying the function to every numerical column(2)
#
write.csv(shapiro_result,file='shapiro_test_result.csv',
          row.names = colnames(cui[-c(1)]))#exporting the output
#
print(shapiro_result)

plot(corr_coef(cui,method = c('spearman')))
corr



