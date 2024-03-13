# Script for testing for differential abundnace  between sexes by using the R package corncob
library(corncob)
library(tidyverse)

set.seed(100)


# Prepare data
# Use tha phyloseq obejct and define which comparisons you want to test as follows
# subset the object for the spesific comparison and for Family level

# NOTE: before that make sure that sam data are factors!! Run this to do convert all meta data to factors
df <- as.data.frame(lapply(sample_data(data.full),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(data.full)
sample_data(data.full) <- sample_data(df)

# create the comparison obj
# column_of_interest is the name of the column in the metada // sam_data

Comparison_1 <- data.full %>% phyloseq::subset_samples(sex %in% c("F","M")) %>%
  tax_glom("Family") 

# Run the test
 Results_comparsion_1red <-differentialTest(formula = ~ sex,
                                          phi.formula = ~ sex,
                                          formula_null = ~ 1,
                                          phi.formula_null = ~ sex,
                                          test = "Wald", boot = FALSE, # LRT is the alternative test you could use 
                                          data = Comparison_1red,
                                          fdr_cutoff = 0.05)

  # If you want to jointly test for differential abundance and differential variability across the 'column of interest'
  # Results_comparsion_1red_b <-differentialTest(formula = ~ column_of_interest,
  #                                           phi.formula = ~ column_of_interest,
  #                                           formula_null = ~ 1,
  #                                           phi.formula_null = ~ 1,
  #                                           test = "Wald", boot = FALSE, # LRT
  #                                           data = Comparison_1,
  #                                           fdr_cutoff = 0.05)


# Save statistics as tables
# significant otus
otuCorn <-otu_to_taxonomy(OTU = Results_comparsion_1red$significant_taxa, data = Comparison_1)
otuCorn<- as.data.frame(otuCorn)
otuCorn<-tibble::rownames_to_column(otuCorn,"OTU.x")

# p values and p adj values
pCorncob <-as.data.frame(Results_comparsion_1red$p,rownames=NULL)
pCorncob<-tibble::rownames_to_column(pCorncob,"OTUp")
colnames(pCorncob)[2]<- c("p")

# p adj
padjCorncob <-as.data.frame(Results_comparsion_1red$p_fdr,rownames=NULL)
padjCorncob<-tibble::rownames_to_column(padjCorncob,"OTUpadj")
colnames(padjCorncob)[2]<- c("padj")

# merge all together 
stats_comparison_1 <-merge(pCorncob,otuCorn, by.x="OTUp",by.y="OTU.x")
stats_comparison_1 <-merge(padjCorncob,stats_comparison_1, by.x="OTUpadj",by.y="OTUp")
stats_comparison_1[order(stats_comparison_1[,2], decreasing = F), ] # ordert them

# save stat
write.csv(stats_comparison_1,"all_stats_corncob_ReducedDataset_sex.csv")

# significant models
Comparison_1_all_models <- Results_comparsion_1red$significant_models


# save them 
cat(capture.output(print(Comparison_1_all_models), file="corncob_ReducedDataset_sex_all_models.txt"))
