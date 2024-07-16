# load libraries
library(ggplot2)
library(hrbrthemes)

# Copy cleaned_trees.rds from tree_and_traits repo
# Copy species_with_ranges.rds from gbif-bulk repo

# Read in trees
cleaned_trees <- readRDS("data/raw-data/cleaned_trees.rds")

# Read in ranges
ranges <- readRDS("data/raw-data/species_with_ranges.rds")

####
# ---- Range sizes per family ----
####

pdf("figures/range_size_distributions.pdf")

for(i in 1:length(unique(ranges$gbif_family))){
  
  tmp <- ranges[ranges$gbif_family==unique(ranges$gbif_family)[i],]
  
  print(
  ggplot(tmp, aes(x=log(aoo))) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + 
    xlim(0,max(log(ranges$aoo))+1) +
    ggtitle(paste(unique(ranges$gbif_family)[i],"| n =",length(tmp[,1])))
  )
  
  if(i == 1){
    
    range_stats_fam <- data.frame(unique(tmp$gbif_family),mean(tmp$aoo), median(tmp$aoo),sd(tmp$aoo),length(na.omit(tmp$aoo)))
    
    colnames(range_stats_fam)<-c("family","mean","median","sd","n")
    
  } else {
    
    range_stats_tmp <- data.frame(unique(tmp$gbif_family),mean(tmp$aoo), median(tmp$aoo),sd(tmp$aoo),length(na.omit(tmp$aoo)))
    
    colnames(range_stats_tmp)<-c("family","mean","median","sd","n")
    
    range_stats_fam <- rbind(range_stats_fam,range_stats_tmp)
    
  }
  
  
}

dev.off()

hist(range_stats_fam$mean)

####
# ---- Ridgeplot ----
####

ranges_aoo <- ranges[!is.na(ranges$aoo),]
ranges_aoo <- droplevels(ranges_aoo)


# library
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(dplyr)
library(forcats)

# Plot
ranges_aoo %>%
  mutate(gbif_family = fct_reorder(gbif_family, aoo, .fun='median')) %>%
  ggplot(aes(x = log(aoo), y = gbif_family, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  labs(title = 'Range sizes per family') +
  theme_minimal() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )

ggsave("figures/ridgeplot_log_AOO_per_family.pdf",width=15,height=30)

####
# ---- Plot range sizes on phylogeny ----
####

####
# ---- Plot range size on phylogenetic tree
####

library(ggtree)
library(TDbook)


#add underscores to match phylo tip labels
aoo<-ranges[,c("species_name","aoo")]
aoo$species_name<-gsub(" ","_",aoo$species_name)
aoo$aoo<-log(aoo$aoo)
colnames(aoo)[2]<-"log_aoo"

pdf("figures/range_size_on_trees.pdf")

for(i in 1:length(cleaned_trees)){
  
  tree_range <- aoo[aoo$species_name%in%cleaned_trees[[i]]$tip.label,]

    #make ggtree
  p <- ggtree(cleaned_trees[[i]])
  #p + geom_tiplab(align=TRUE,size=1)
  
  ## visualize  Trait data using bar chart
  ## and align based on tree structure
  print(
  p + geom_facet(panel = "Range", data = tree_range , geom = geom_col, 
                 aes(x = log_aoo), orientation = 'y', width = .6) +
    theme_tree2(legend.position=c(.05, .85))
  )
  
}
  
dev.off()
