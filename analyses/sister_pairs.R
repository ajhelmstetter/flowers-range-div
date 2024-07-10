# load libraries
library(ggplot2)
library(hrbrthemes)
library(phytools)
library(phangorn)

# Copy cleaned_trees.rds from tree_and_traits repo
# Copy species_with_ranges.rds from gbif-bulk repo

# Read in trees
cleaned_trees <- readRDS("data/raw-data/cleaned_trees.rds")

# Read in ranges
ranges <- readRDS("data/raw-data/species_with_ranges.rds")

pdf("figures/sister_pair_age_vs_range_size_diff.pdf")

for(i in 1:length(cleaned_trees)){
  
  # remove underscores so names match
  cleaned_trees[[i]]$tip.label <- gsub("_"," ", cleaned_trees[[i]]$tip.label)
  
  #make df of species with range data
  tmp_ranges <- ranges[ranges$species_name %in% cleaned_trees[[i]]$tip.label,]
  tmp_ranges <- tmp_ranges[,c("species_name", "range_size")]
  
  # get species without ranges
  na_df <- data.frame(setdiff(cleaned_trees[[i]]$tip.label, tmp_ranges$species_name), rep(NA, length(setdiff(cleaned_trees[[i]]$tip.label, tmp_ranges$species_name))))
  colnames(na_df)<-colnames(tmp_ranges)
  
  #combined dfs
  #tmp_ranges <- merge(tmp_ranges, na_df, by ='species_name')
  tmp_ranges <- tmp_ranges[match(cleaned_trees[[i]]$tip.label,tmp_ranges$species_name),]
  table(tmp_ranges$species_name,cleaned_trees[[i]]$tip.label)
  
  phy <- cleaned_trees[[i]]
  
  # get all sister pairs
  is.ultrametric(phy)
  sisters <- diverge::extract_sisters(tree=phy, sis_age=TRUE)
  
  sisters$range_sp1<-rep(NA,length(sisters[,1]))
  sisters$range_sp2<-rep(NA,length(sisters[,1]))
  sisters$range_size_diff<-rep(NA,length(sisters[,1]))
  
  for(j in 1:length(sisters[,1])){
    
    #get range sizes for each of sisters
    
    if(length(tmp_ranges$range_size[grep(sisters[j,1], x = tmp_ranges$species_name)]) == 1){
      
      sisters[j,"range_sp1"] <- tmp_ranges$range_size[grep(sisters[j,1], x = tmp_ranges$species_name)]
      
    } else {
      
      sisters[j,"range_sp1"] <- NA
      
    }
    
    if(length(tmp_ranges$range_size[grep(sisters[j,2], x = tmp_ranges$species_name)]) == 1){
      
      sisters[j,"range_sp2"] <- tmp_ranges$range_size[grep(sisters[j,2], x = tmp_ranges$species_name)]
      
    } else {
      
      sisters[j,"range_sp2"] <- NA
      
    }
    
    
    # get divergence times
    
    sisters$range_size_diff[j] <- max(sisters$range_sp1[j],sisters$range_sp2[j])-min(sisters$range_sp1[j],sisters$range_sp2[j])
    
  }
  
  sisters
  
  # make linear model
  #mod <- summary(lm(sisters$pair_age~log(sisters$range_size_diff)))
  
  #NOTE: not sure when to log range diff
  print(
  ggplot(sisters, aes(x=pair_age, y=log(range_size_diff))) +
    geom_point(size=1,shape=21,alpha=0.5) +
    geom_smooth(method=lm , color="purple", alpha=0.4, fill="lavender", se=TRUE) + 
    annotate("text", x = -Inf, y = Inf, label = names(cleaned_trees)[i], hjust = -.2, vjust = 2) +
    annotate("text", x = -Inf, y = Inf, label = paste("Sister pairs =", sum(!is.na(sisters$range_size_diff))), hjust = -.2, vjust = 4) +
    #annotate("text", x = -Inf, y = Inf, label = paste("adj r2 =", round(mod$adj.r.squared,3),"p-value =",round(mod$coefficients[2,4],3)), hjust = -.2, vjust = 6) +
    theme_minimal()
  )
  
  
  if(i == 1){
    
    combined_df<-sisters
    
  } else {
    
    combined_df <- rbind(combined_df, sisters)
    
  }
  

}

dev.off()

#remove trees with extreme DR values
#combined_df<-combined_df[combined_df$tree!="FabalesRosalesFagalesCucurbitales-Kates_et_al-2024.tree",]
#combined_df<-combined_df[combined_df$tree!="Dipterocarpaceae-Bansal_et_al-2022.tree",]

# make linear model
mod <- summary(lm(combined_df$pair_age~log(combined_df$range_size_diff)))

ggplot(combined_df, aes(x=pair_age, y=log(range_size_diff))) +
  geom_point(size=1,shape=21,alpha=0.5) +
  geom_smooth(method=lm , color="purple", alpha=0.4, fill="lavender", se=TRUE) + 
  annotate("text", x = -Inf, y = Inf, label = paste("No. sisters = ", length(combined_df[,1])), hjust = -.4, vjust = 4) +
  annotate("text", x = -Inf, y = Inf, label = paste("adj r2 =", round(mod$adj.r.squared,3),"p-value =",round(mod$coefficients[2,4],3)), hjust = -.2, vjust = 6) +
  theme_minimal()

