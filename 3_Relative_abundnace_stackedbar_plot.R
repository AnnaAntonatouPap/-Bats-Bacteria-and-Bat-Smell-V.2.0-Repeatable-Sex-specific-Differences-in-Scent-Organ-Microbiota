# Create stacked bar plots with relative abundance in the Family level, showing in different colors only the Top 20 families and the rest pooled together 

# stacked bar plot with all the samples 

prune.dat <- prune_taxa(taxa_sums(data.full) > 2, data.full.filt) # remove less than 2%
top20 <- names(sort(taxa_sums(prune.dat), decreasing=TRUE)[1:43])
#rest <- names(sort(taxa_sums(prune.dat), decreasing=TRUE)[20:277])

taxtab20 = cbind(tax_table(prune.dat), Family_top20 = NA)
taxtab20[top20, "Family_top20"] <- as(tax_table(prune.dat)[top20, "Family"],
                                        "character")
tax_table(prune.dat) <- tax_table(taxtab20)

prune.dat <-transform_sample_counts(prune.dat, function(x) x/sum(x))
data.frame = psmelt(prune.dat)

# colors
stackColors1 <- c("#673770", "#5F7FC7", "orange","#DA5724", "#5e738f", "#CD9BCD",
                 "#AD6F3B", "#CBD588","#D14285", "#508578" , "#C84248", "#8569D5", 
                 "#652926","#D1A33D", "bisque3", "#599861","red","blue","pink","cyan2")


# Plot
plot<-ggplot(data.frame, aes(x=Sample, y=Abundance, fill=Family_top20)) + 
  geom_bar(stat="identity") + # position = "fill"
  theme_bw() +
  scale_fill_manual(values=stackColors1)+
  facet_grid(~sex, scale="free")+ #My_Theme +
  theme(text=element_text(size=20,family="Times New Roman")) +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),axis.title.x = element_blank(), axis.title.y =       element_text(size=20), axis.text.y = element_text(size=13)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  ylab("Relative abundance") +labs(fill="Family")+theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 13,face = "italic"),      # Adjust the text size in the legend
  legend.title = element_text(size = 13))
plot


# stacked bar for Males 
data_males <-subset_samples(data.full, sex=="M")

prune.dat_males <- prune_taxa(taxa_sums(data_males) > 2, data_males) #remove less than 2%
top20males <- names(sort(taxa_sums(prune.dat_males), decreasing=TRUE)[1:48])

taxtab20 = cbind(tax_table(prune.dat_males), Family_top20 = NA)
taxtab20[top20males, "Family_top20"] <- as(tax_table(prune.dat_males)[top20males, "Family"],
                                        "character")
tax_table(prune.dat_males) <- tax_table(taxtab20)

prune.dat_males <-transform_sample_counts(prune.dat_males, function(x) x/sum(x))
data.frame_males = psmelt(prune.dat_males)

# colors
stackColors1 <- c("#673770", "#5F7FC7", "orange","#DA5724", "#5e738f", "#CD9BCD",
                 "#AD6F3B", "#CBD588","#D14285", "#508578" , "#C84248", "#8569D5", 
                 "#652926","#D1A33D", "bisque3", "#599861","red","blue","pink","cyan2"
)


plot<-ggplot(data.frame_males, aes(x=Sample, y=Abundance, fill=Family_top20)) + 
  geom_bar(stat="identity") + # position = "fill"
  theme_bw() +
  scale_fill_manual(values=stackColors1)+
  #facet_grid(~Sample_type+sampling_period, scale="free")+
  theme(text=element_text(size=12,family="Times New Roman"))+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),axis.title.x = element_blank(), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  ylab("Relative abundance") +labs(fill="Family")+theme(legend.position = "bottom")
plot

# stacked bar plot females
data_females <-subset_samples(data.full, sex=="F")

prune.dat_females <- prune_taxa(taxa_sums(data_females) > 2, data_males) #remove less than 2%
top20males <- names(sort(taxa_sums(prune.dat_females), decreasing=TRUE)[1:48])

taxtab20 = cbind(tax_table(prune.dat_females), Family_top20 = NA)
taxtab20[top20males, "Family_top20"] <- as(tax_table(prune.dat_males)[top20males, "Family"],
                                        "character")
tax_table(prune.dat_females) <- tax_table(taxtab20)

prune.dat_females <-transform_sample_counts(prune.dat_females, function(x) x/sum(x))
data.frame_females = psmelt(prune.dat_females)

# colors
stackColors1 <- c("#673770", "#5F7FC7", "orange","#DA5724", "#5e738f", "#CD9BCD",
                 "#AD6F3B", "#CBD588","#D14285", "#508578" , "#C84248", "#8569D5", 
                 "#652926","#D1A33D", "bisque3", "#599861","red","blue","pink","cyan2")


plot<-ggplot(data.frame_females, aes(x=Sample, y=Abundance, fill=Family_top20)) + 
  geom_bar(stat="identity") + # position = "fill"
  theme_bw() +
  scale_fill_manual(values=stackColors1)+
  #facet_grid(~Sample_type+sampling_period, scale="free")+
  theme(text=element_text(size=12,family="Times New Roman"))+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),axis.title.x = element_blank(), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  ylab("Relative abundance") +labs(fill="Family")+theme(legend.position = "bottom")
plot



