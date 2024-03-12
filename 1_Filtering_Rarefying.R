# Importing files, filtering samples, rarafy data, CSS data Bat project
# Load libraries 
library(xml2)
library(XML,lib.loc = "C:/Program Files/R/R-4.0.2/library")
library(rlang, lib.loc = "C:/Program Files/R/R-4.0.2/library")
library(ggplot2)
library(vegan)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(ggpval)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(biomformat)
library(microbiome)
library(DESeq2)
library(dunn.test)
library(rstatix)
library(ggpubr, lib.loc="~/R/win-library/4.0")
library(metagenomeSeq)
library(ape)
library(rfm)
library(picante)
library("dbplyr", lib.loc="~/R/win-library/4.0")
library("geosphere", lib.loc="~/R/win-library/4.0")
library(ggpubr)
library(moments)
library(lme4)
library(car)
library("sjPlot", lib.loc="~/R/win-library/4.0")
library("lme4", lib.loc="~/R/win-library/4.0")
library(broom)
library("multcomp", lib.loc="~/R/win-library/4.0")
library(ggeffects)
library("MuMIn", lib.loc="~/R/win-library/4.0")
library(ecodist)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(tidyverse)
library(extrafont)
font_import()
loadfonts(device = "win") 

# My theme for plotting 
My_Theme = theme(
  axis.title.x = element_text(size = 22,family="Times New Roman",face = "bold", colour = "black"),
  axis.text.x = element_text(size = 20,family="Times New Roman"),
  axis.text.y = element_text(size = 20,family="Times New Roman"),
  axis.title.y = element_text(size = 22,family="Times New Roman",face = "bold", colour = "black"),
  )

# Loading raw files (biom file, meta data as csv and tree as nwk), after this we create the phyloseq object with the name "data"

# Filter out the samples labelled as “NO” (data do not belong to this project) in “new full_dataset” column. The number of remaining samples should be 54.
  this is the dataset we gonna work on, called data.full

######################################################################################################################################################

# Rename Ranks 
colnames(tax_table(data.full)) <- c("Kingdom","Phylum","Class",
                                     "Order","Family","Genus")

# Discard taxa that are classified as Mitochondria or/and Chloroplasts 
# NOTE : no Archea found in this dataset

data.full <- data.full %>%
  subset_taxa(
    ((Family != "Mitochondria") | is.na(Family)) & (Order  != "Chloroplast")


######################################################################################################################################################   
# Rarefy full dataset

set.seed(1)
data.full.rare <- rarefy_even_depth(data.full,rngseed=1, sample.size=min(sample_sums(data.full)), replace=F)
ntaxa(data.full.rare) ## 70 














