
# log transform the phylosew object prior to Beta diversity analysis

 data.full.log <- microbiome::transform(data.full, 'log10p')


# Calculate Beta diversity distance matrices 

# bRAY Curtis
Bray_pcoa <- ordinate(
  physeq = data.reduced.css, 
  method = "PCoA", 
  distance = "bray")
# Ordinate NMDS bray 0.1953305  
Bray_nmds  <- ordinate(
  physeq = data.reduced.css, 
  method = "NMDS", 
  distance = "bray")
# ordinate jaccard NMDS 0.2430796  

# Jaccard
Jaccard_nmds <- ordinate(
  physeq = data.reduced.css, 
  method = "NMDS", 
  distance = "jaccard",binary= TRUE)
# ordinate Jaccard PCoA
Jaccard_pcoa <- ordinate(
  physeq = data.reduced.css, 
  method = "PCoA", 
  distance = "jaccard", binary=TRUE)

# Unweighted unifracs
# ordinate Unifrac PCoA
unifrac_pcoa <- ordinate(
  physeq = data.reduced.css, 
  method = "PCoA", 
  distance = "unifrac")
# ordinate Unifrac NMDS 0.1714688  
unifrac_nmds <- ordinate(
  physeq = data.reduced.css, 
  method = "NMDS", 
  distance = "unifrac")

# Weigheted unifrac  
# ordinate wUnifrac PCoA
wunifrac_pcoa<- ordinate(
  physeq = data.reduced.css, 
  method = "PCoA", 
  distance = "wunifrac")
# ordinate wUnifrac NMDS 0.1670089 
wunifrac_nmds  <- ordinate(
  physeq = data.reduced.css, 
  method = "NMDS", 
  distance = "wunifrac")


  ############################################################################################################################################################################

 # Ordination plot for each distance matrix as follows

plot <-plot_ordination(
  physeq = data.full.log,
  ordination = Bray_pcoa,
  color = "sex",
  title = "") + 
  geom_point (size=5)+
  scale_color_manual(values = c("#660099", "chartreuse4"),name="", breaks=c("F","M"),
                     labels=c("Female","Male"))+
  
  theme_bw() + My_Theme +
  ggtitle("Bray-Curtis")+
  theme(legend.background = element_rect(size=0.3,linetype="solid", colour ="black"))+                      
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot

# same as above for all the plots

##############################################################################################################################################################################33####

# Prepare data for Permanova
# Calculate distance matrices and save them if needed 

otus <-as(otu_table(data.full.log),"matrix")  # otu table as matrix
otus <-t(otus) # transpose
otus <-as.data.frame(otus)

# bray curtis
dist.bray.full = vegdist(otus, method = "bray")
dist.bray.full.matrix <- as.matrix(dist.bray.full)
#write.csv(dist.bray.full.matrix,"BrayCurtis.csv")

# jaccard
dist.jaccard.full = vegdist(otus.css, method = "jaccard",binary=T)
dist.jaccard.full.matrix <- as.matrix(dist.jaccard.full)
write.csv(dist.jaccard.full.matrix,"Jaccard_fullDataset_log.csv")

# weighted unifrac
dist.Wunifrac.full <- UniFrac(data.full.log, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)
dist.Wunifrac.full.matrix <- as.matrix(dist.Wunifrac.full)

# unweighted unifrac 
dist.UNunifrac.full <- UniFrac(data.full.log, 
                                        weighted = FALSE, 
                                        normalized = FALSE,  
                                        parallel = FALSE, 
                                        fast = TRUE)

dist.UNunifrac.full.matrix <- as.matrix(dist.UNunifrac.full)

# Betadipser test for sex,  for Bray Curtis ( do the same for the rest distnace matrices  and for area too)

betadisper_Bray <- betadisper(dist.bray.full,data.full.log@sam_data$sex)

# ANOVA or PERMTEST
anova(betadisper_Bray)
perm_betadisper_Bray <- permutest(betadisper_Bray, permutations = 9999, pairwise = TRUE)


# PERMANOVA example for Bray Curtis// Run similar for the rest

permanova_bray_curtis<-adonis2(dist.bray.full.matrix ~ sex + area,
                            by=NULL, 
                            permutations = 9999,
                             method = "bray",
                             data = data.full.log@sam_data)
