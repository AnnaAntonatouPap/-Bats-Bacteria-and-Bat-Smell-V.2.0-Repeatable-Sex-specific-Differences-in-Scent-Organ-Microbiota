# Importing files, filtering samples, rarafy data, CSS data Bat project
# Load libraries 
library(xml2)
library(XML,lib.loc = "C:/Program Files/R/R-4.0.2/library")
library(rlang, lib.loc = "C:/Program Files/R/R-4.0.2/library")
library(ggplot2)
library(vegan)
library(dplyr)
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




# Loading raw files (biom file, meta data as csv and tree as nwk), after this we create the phyloseq object with the name 
# data full



# Rename Ranks 

colnames(tax_table(data.full)) <- c("Kingdom","Phylum","Class",
                                     "Order","Family","Genus")

# Discard taxa that are classified as Mitochondria or/and Chloroplasts 
# NOTE : no Archea found in this dataset

data.full.filt <- data.full %>%
  subset_taxa(
    ((Family != "Mitochondria") | is.na(Family)) & (Order  != "Chloroplast")
  )





# css norm fulldataset


library(metagenomeSeq)

## Convert the phyloseq object to a metagenomeSeq object (MRexperiment).
metagenome.obj <- phyloseq_to_metagenomeSeq(data.full.filt)
## Calculate the proper percentile by which to normalize counts
cNstat <- metagenomeSeq::cumNormStatFast(metagenome.obj) 
# cNstat  #0.5
## Normalise counts
metagenome.obj <- metagenomeSeq::cumNorm(metagenome.obj, p = cNstat)
## Export the normalised count table
metag.norm.counts <- metagenomeSeq::MRcounts(metagenome.obj, norm = TRUE)
## Add a pseudocount of 0.0001 to the table and log transform
metag.norm.counts_log <- log(metag.norm.counts+0.0001)
## Substract the value from log of pseudocount to preserve zeros of the original counts
metag.norm.counts_log2 <- metag.norm.counts_log-(log(0.0001))

## Make a new phyloseq object with with the new OTU table
otu_normMG.obj <- otu_table(metag.norm.counts_log2, taxa_are_rows = TRUE)
phylo_normMG.obj <- phyloseq(otu_normMG.obj, data.full.filt@tax_table, data.full.filt@sam_data)
data.full.css <- merge_phyloseq(phylo_normMG.obj, data.full.filt@phy_tree)

# save OTU
OTUmatrix =as(otu_table(data.full.css),"matrix") 
OTUdf = as.data.frame(OTUmatrix)
write.csv(OTUdf, "FullDataset_css_transformed.csv")

# Rarefy full dataset

set.seed(1)
data.full.rare <- rarefy_even_depth(data.full.filt,rngseed=1, sample.size=min(sample_sums(data.full.filt)), replace=F)
ntaxa(data.full.rare) ## 70 

# Rarefy reduced dataset

set.seed(1)
data.reduced.rare <- rarefy_even_depth(data.reduced.filt,rngseed=1, sample.size=min(sample_sums(data.reduced.filt)), replace=F)
ntaxa(data.reduced.rare) ## 114 


# violin plot for read count

# extract otu table of the filtered full dataset by sex

otu_violin = as(otu_table(data.full.filt), "matrix")
if(taxa_are_rows(data.full.filt)){otu_violin <- t(otu_violin)} # transpose

otu_violindf = as.data.frame(otu_violin) # make dataframe

# total read count per sample
otu_violindf <-otu_violindf %>%
     mutate(Total = rowSums(.))

# keep total and samples
samples_total  <- otu_violindf %>% dplyr::select("Total")

# paste metadata sex
meta_violin = as(sample_data(data.full.filt), "matrix")
meta_violin = as.data.frame(meta_violin)

meta_violin <- meta_violin %>% dplyr::select("sex")

samples_total <-cbind(samples_total, sex= meta_violin$sex)

# violin plot prep
samples_total$sex <- as.factor(samples_total$sex)
# sample size
sample_size = samples_total %>% group_by(sex) %>% summarize(num=n())

# save 
write.csv(samples_total,"Total_counts_per_sample_fullDataset.csv")


# Violin plot for counts 

dp <- ggplot(samples_total, aes(x=sex, y=Total, fill=sex)) + 
  geom_violin(trim=T,alpha=0.7,adjust=2,show.legend = FALSE)+
  geom_boxplot(width=0.05, fill="white", show.legend = FALSE)+
  geom_point(position = position_jitter(seed = 1, width = 0.2), show.legend = FALSE) +
  #My_Theme+
  theme_bw() + My_Theme +
  #theme(text=element_text(size=12,family="Times New Roman"))+
  #theme(axis.text.x = element_blank(),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12))+
   xlab("Sex") + ylab("Total counts per sample")

dp<-dp+scale_fill_manual(values=c("#660099", "chartreuse4"))#+ theme(legend.title=element_blank())
dp

ggsave("violinplot_reads_fullDataset.pdf",dp,device = cairo_pdf,width=400,height = 350, units ="mm",dpi=1000)

# plot boxplot
dp2 <- ggplot(samples_total, aes(x=sex, y=Total, fill=sex)) +
    geom_boxplot() +
    geom_point(position = position_jitter(seed = 1, width = 0.2), show.legend = FALSE) +
    #geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_bw() + My_Theme +
   xlab("Sex") + ylab("Total counts per sample")

dp2<-dp2+scale_fill_manual(values=c("#ae5a41", "#006666"))#+ theme(legend.title=element_blank())
dp2



#ggsave("boxplot_reads_fullDataset.pdf",dp2,device = cairo_pdf,width=400,height = 350, units ="mm",dpi=1000)

# add stats
library(ggpubr)

dp2+stat_compare_means()
dp2

# run the test
compare_means(Total ~ sex, data = samples_total)


ggsave("boxplot_reads_fullDataset_withtest_no_dots.pdf",last_plot(),device = cairo_pdf,width=400,height = 350, units ="mm",dpi=1000)



# recurve for fulldataset
###3 rarecurve by sample for full dataset
cairo_pdf("accumulation_curve_full.pdf", fallback_resolution=300 ) # done

rarecurve(t(otu_table(data.full.filt)), cex=0.5)
dev.off()


