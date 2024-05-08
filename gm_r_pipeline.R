#WRKY CUI analysis

#importing necessary libraries
library(tidyverse)
library(metan)
library(corrplot)
library(readxl)
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
cui_df <- joined_df[-c(7:8,12:18)]
view(cui_df)


