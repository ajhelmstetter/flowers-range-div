# Copy cleaned_trees.rds from tree_and_traits repo

# Read in trees
cleaned_trees <- readRDS("data/raw-data/cleaned_trees.rds")

# loop through trees and calculate DR for each species

for(i in 1:length(cleaned_trees)){
  
  if(i == 1){
    
    dr <- epm::DRstat(cleaned_trees[[i]])
    
    dr_df <- data.frame(rep(names(cleaned_trees[i]),length(names(dr))), names(dr), dr)
    
    colnames(dr_df) <- c("tree", "species", "dr")
    
  } else {
    
    dr <- epm::DRstat(cleaned_trees[[i]])
    
    tmp <- data.frame(rep(names(cleaned_trees[i]),length(names(dr))), names(dr), dr)
    
    colnames(tmp) <- c("tree", "species", "dr")
    
    dr_df <- rbind(dr_df, tmp)
    
  }
      
}

head(dr_df)
str(dr_df)

#write DRs to file
saveRDS(dr_df,"outputs/dr_df.rds")

