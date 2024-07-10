# load libraries
library(ggplot2)
library(hrbrthemes)

# Copy cleaned_trees.rds from tree_and_traits repo
# Copy species_with_ranges.rds from gbif-bulk repo

# Read in trees
cleaned_trees <- readRDS("data/raw-data/cleaned_trees.rds")

# Read in ranges
ranges <- readRDS("data/raw-data/species_with_ranges.rds")

# Read in diversification rates
dr_df <- readRDS("outputs/dr_df.rds")

#remove '_' from tip labels
dr_df$species <- gsub("_"," ",dr_df$species)

####
# ---- Plot per-tree range size vs DR ----
####

pdf("figures/range_size_vs_dr.pdf")

for(i in 1:length(cleaned_trees)){
  
  tree_dr <- dr_df[dr_df$tree==unique(dr_df$tree)[i],]
  
  tree_dr <- tree_dr[tree_dr$species%in%ranges$species_name,]
  tree_dr <- tree_dr[order(tree_dr$species),]
  
  tree_range <- ranges[ranges$species_name%in%tree_dr$species,]
  tree_range <- tree_range[order(tree_range$species_name),]
  
  # NON-MERGE SOLUTION WITH DUPLICATES ISSUE
  # # Check matches are good
  # print(table(tree_dr$species==tree_range$species_name))
  # 
  # tree_range_df<-data.frame(tree_dr$species,tree_dr$dr,log(tree_range$range_size))
  # colnames(tree_range_df) <- c("species","dr","log_range_size")
  
  #make column names match
  colnames(tree_range)[1] <- "species"
  
  #duplicates will still be included
  tree_range_df <- merge(tree_range,tree_dr,by="species")
  tree_range_df$log_range_size <- log(tree_range_df$range_size)
  
  # make linear model
  mod <- summary(lm(tree_range_df$log_range_size~tree_range_df$dr))
  
  # linear trend + confidence interval
    print(
    ggplot(tree_range_df, aes(x=log_range_size, y=dr)) +
    geom_point(size=1,shape=21,alpha=0.5) +
    geom_smooth(method=lm , color="darkblue", alpha=0.4, fill="#69b3a2", se=TRUE) + 
    annotate("text", x = -Inf, y = Inf, label = unique(dr_df$tree)[i], hjust = -.2, vjust = 2) +
    annotate("text", x = -Inf, y = Inf, label = paste("No. species = ", length(unique(tree_range_df$species))), hjust = -.4, vjust = 4) +
      annotate("text", x = -Inf, y = Inf, label = paste("adj r2 =", round(mod$adj.r.squared,3),"p-value =",round(mod$coefficients[2,4],3)), hjust = -.2, vjust = 6) +
    theme_minimal()
  )

    if(i == 1){
      
      combined_df<-tree_range_df
      
    } else {
      
      combined_df <- rbind(combined_df, tree_range_df)
      
    }
    
    
}

dev.off()


#remove trees with extreme DR values
#combined_df<-combined_df[combined_df$tree!="FabalesRosalesFagalesCucurbitales-Kates_et_al-2024.tree",]
#combined_df<-combined_df[combined_df$tree!="Dipterocarpaceae-Bansal_et_al-2022.tree",]

# make linear model
mod <- summary(lm(combined_df$log_range_size~combined_df$dr))

ggplot(combined_df, aes(x=log_range_size, y=log(dr))) +
  geom_point(size=1,shape=21,alpha=0.5) +
  geom_smooth(method=lm , color="darkblue", alpha=0.4, fill="#69b3a2", se=TRUE) + 
  annotate("text", x = -Inf, y = Inf, label = paste("No. species = ", length(unique(combined_df$species))), hjust = -.4, vjust = 4) +
  annotate("text", x = -Inf, y = Inf, label = paste("adj r2 =", round(mod$adj.r.squared,3),"p-value =",round(mod$coefficients[2,4],3)), hjust = -.2, vjust = 6) +
  theme_minimal()
