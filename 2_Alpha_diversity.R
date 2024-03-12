# Calculate alpha diveristies for unraraefied data

# A. Shannon, Chao and Observed 
alpha.full_unrare <- estimate_richness(data.full, 
                                        measures = c("Shannon","Chao1","Observed"))

# B. Faith's PD
# Prepare data 
pd_otu <- as.data.frame(data.full.filt@otu_table) #otu table as data frame
pd_tree <- data.full.filt@phy_tree # tree
pd_tree=multi2di(pd_tree) # OPTIONAL, if tree is NOT rooted,  convert it to rooted

# Calculate phylogenetic diversity
df.pd <- pd(t(pd_otu),pd_tree,include.root=T) 


# All diversities in one dataframe with all metadata together
all_alpha_full_unrare <- cbind(sample_data(data.full.filt), alpha.full_unrare,df.pd)


# In order to run LMs or LMMs we need to check if the data are normally distributed
# To do so we use the Shapiro test and we check normality for all three diversities


shapiro.test(all_alpha_full_unrare$Shannon) #NOT 0.001452
shapiro.test(all_alpha_full_unrare$Chao1)# NOT  0.0007271
shapiro.test(all_alpha_full_unrare$PD)
shapiro.test(all_alpha_full_unrare$Observed)

# If any of them is not normaly distributed, we transform tha data to reach nromality 

# Square root transformation of  Oberved and Chao1 and check of normality  
shapiro.test(sqrt(all_alpha_full_unrare$Observed)) 
shapiro.test(sqrt(all_alpha_full_unrare$Chao1)) 

# Include the new transformed values of Observed and Chao1 in the  big dataframe 
all_alpha_full_unrare$Observed.trans <- sqrt(all_alpha_full_unrare$Observed)
all_alpha_full_unrare$Chao.trans <- sqrt(all_alpha_full_unrare$Chao1)



# Plot diversities// Boxplots for sex 


my_comparisons = list( c("F", "M"))

# Plot Shannon
dp2 <- ggplot(all_alpha_full_unrare, aes(x=sex, y=Shannon, fill=sex)) +
    geom_boxplot() +
    theme_bw() + My_Theme  + guides(fill = 'none') +
   xlab("Sex") + ylab("Shannon's diversity")

dp2<-dp2+scale_fill_manual(values=c("#ff8800", "#004fff"))

dp2 + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label="p.signif")


# Plot Observed No.of OTUs

dp3 <- ggplot(all_alpha_full_unrare, aes(x=sex, y=Observed.trans, fill=sex)) +
    geom_boxplot() +
    theme_bw() + My_Theme  + guides(fill = 'none') +
   xlab("Sex") + ylab("Observed No of OTUs")

dp3<-dp3+scale_fill_manual(values=c("#ff8800", "#004fff"))
dp3 +stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label="p.signif")

# Plot Chao1

dp4 <- ggplot(all_alpha_full_unrare, aes(x=sex, y=Chao.trans, fill=sex)) +
    geom_boxplot() +
    theme_bw() + My_Theme + guides(fill = 'none') +
   xlab("Sex") + ylab("Chao1")

dp4<-dp4+scale_fill_manual(values=c("#ff8800", "#004fff")
dp4 +stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label="p.signif")



############################################################################################################################################################################

#  Use Linear Models to investigate  weather alpha diversity metrics are influenced by sex and area 


## Shannon
model1 <- lm(Shannon ~ sex + area, data=all_alpha_full_unrare)
#model1T <- lm(Shannon ~ area + sex, data=all_alpha_full_unrare) # run the same model with different order of the  explanatory variables, no change observed 

# summerize and save as tavle 
summary(model1)
tab_model(model1,show.se=TRUE,string.se = "std. Error")

## Observed 
model2 <- lm(Observed.trans ~ sex + area, data=all_alpha_full_unrare)
#model2T <- lm(Observed.trans ~ area + sex, data=all_alpha_full_unrare) # run the same model with different order of the  explanatory variables, no change observed 

summary(model2)
tab_model(model2,show.se=TRUE,string.se = "std. Error")

## Chao1 
model3 <- lm(Chao.trans ~ sex + area, data=all_alpha_full_unrare)
#model3T <- lm(Chao.trans ~ area + sex, data=all_alpha_full_unrare) # run the same model with different order of the  explanatory variables, no change observed 

# summerize and save as tavle 
summary(model3)
tab_model(model3,show.se=TRUE,string.se = "std. Error")

## faith's PD 

model4 <- lm(PD ~ sex + area, data=all_alpha_full_unrare)
#model4T <- lm(PD.trans ~ area + sex, data=all_alpha_full_unrare) # run the same model with different order of the  explanatory variables, no change observed 

# summerize and save as tavle 

summary(model4)
tab_model(model4,show.se=TRUE,string.se = "std. Error")
