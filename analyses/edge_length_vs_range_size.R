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

pdf("figures/terminal_branch_length_vs_range_size.pdf")

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
  
  # http://blog.phytools.org/2013/10/finding-edge-lengths-of-all-terminal.html
  n<-length(phy$tip.label)
  edge_length<-setNames(phy$edge.length[sapply(1:n,function(x,y)   which(y==x),y=phy$edge[,2])],phy$tip.label)
  edge_length
  
  tmp_ranges$species_name == names(edge_length)

  ranges_bl <- cbind(tmp_ranges,edge_length)
  
  # make linear model
  mod <- summary(lm(log(ranges_bl$range_size)~ranges_bl$edge_length))
  
  #NOTE: not sure when to log range diff
  print(
    ggplot(ranges_bl, aes(x=edge_length, y=log(range_size))) +
      geom_point(size=1,shape=21,alpha=0.5) +
      geom_smooth(method=lm , color="darkgreen", alpha=0.4, fill="lightgreen", se=TRUE) + 
      annotate("text", x = -Inf, y = Inf, label = names(cleaned_trees)[i], hjust = -.2, vjust = 2) +
      #annotate("text", x = -Inf, y = Inf, label = paste("Sister pairs =", sum(!is.na(sisters$range_size_diff))), hjust = -.2, vjust = 4) +
      annotate("text", x = -Inf, y = Inf, label = paste("adj r2 =", round(mod$adj.r.squared,3),"p-value =",round(mod$coefficients[2,4],3)), hjust = -.2, vjust = 6) +
      theme_minimal()
  )
  
  if(i == 1){
    
    combined_df<-ranges_bl
    
  } else {
    
    combined_df <- rbind(combined_df, ranges_bl)
    
  }
  
}

dev.off()

#remove trees with extreme DR values
#combined_df<-combined_df[combined_df$tree!="FabalesRosalesFagalesCucurbitales-Kates_et_al-2024.tree",]
#combined_df<-combined_df[combined_df$tree!="Dipterocarpaceae-Bansal_et_al-2022.tree",]

# make linear model
mod <- summary(lm(combined_df$pair_age~log(combined_df$range_size_diff)))

ggplot(combined_df, aes(x=edge_length, y=log(range_size))) +
  geom_point(size=1,shape=21,alpha=0.5) +
  geom_smooth(method=lm , color="darkgreen", alpha=0.4, fill="lightgreen", se=TRUE) + 
  annotate("text", x = -Inf, y = Inf, label = paste("No. species  = ", length(unique(combined_df$species_name))), hjust = -.4, vjust = 4) +
  annotate("text", x = -Inf, y = Inf, label = paste("adj r2 =", round(mod$adj.r.squared,3),"p-value =",round(mod$coefficients[2,4],3)), hjust = -.2, vjust = 6) +
  theme_minimal()
