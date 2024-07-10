rm(list=ls())

#load libraries
library(RPANDA)

#running Julia in R
#https://cran.r-project.org/web/packages/JuliaCall/readme/README.html
library(JuliaCall)

####
# ---- ClaDS with Julia in R ----
####

#Issue with Julia in R
#https://github.com/Non-Contradiction/JuliaCall/issues/224
#fixed by adding R_LD_LIBRARY_PATH manually to ~/.zshrc
#add this to end of ~/.zshrc :
#export R_LD_LIBRARY_PATH="/home/andrew.helmstetter/Programs/julia-1.10.2/lib/julia:/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu"
#Then need to open terminal and run
#R -q          
#library("JuliaCall")
#julia <- julia_setup(JULIA_HOME = "/home/andrew.helmstetter/Programs/julia-1.10.2/bin/")
#Could probably edit PATH in R somehow to fix this

#set up julia
julia <- julia_setup(JULIA_HOME = "/home/andrew.helmstetter/Programs/julia-1.10.2/bin/")

#load library
julia_library("PANDA")
julia_library("JLD2")

# Copy cleaned_trees from trees_and_traits repo
system("cp -r ~/Dropbox/projects/AJH_FloweRS/trees_and_traits/data/derived-data/cleaned_trees data/raw-data/")

# clads needs suffix '.tre' to work
system("rename \"s/.tree$/.tre/\" ~/Dropbox/projects/AJH_FloweRS/flower-range-div/data/raw-data/cleaned_trees/*.tre")

#get full file paths to trees
list_paths <- paste(here::here(),"/",list.files("data/raw-data/cleaned_trees", full.names = TRUE),sep="")
list_names <- gsub(".tree","",list.files("data/raw-data/cleaned_trees"))

#read in sampling fractions



#loop through trees
for(i in 1:length(list_paths)){
  
  paste("my_tree = load_tree(\"",list_paths[1],"\")",sep="")
  
  #read in tree
  julia_command("my_tree = load_tree(\"/home/andrew.helmstetter/Dropbox/projects/AJH_FloweRS/flowers-range-div/data/raw-data/cleaned_trees/Allium-Han_et_al-2019.tre\")")

  #run ClaDS
  julia_command("output = infer_ClaDS(my_tree, f = 0.72, print_state = 100)")
  
  #save for manipulation in R
  julia_command("save_ClaDS_in_R(output, \"/home/andrew.helmstetter/Dropbox/projects/AJH_FloweRS/phlox/outputs/clads_out.Rdata\")")
  
  
}

####
# ---- Analyse ClaDS output in R
####

#load ClaDS output
load("outputs/clads_out.Rdata")

#look at output structure
#NOT RUN: some julia error if julia not called
#str(CladsOutput)

#read in tree
#NOTE: ladderize can mess up rate plotting so ladderize before running ClaDS 
#NOTE: (not sure if this actually helps?)
phy <- read.tree("data/derived-data/polemoniaceae_corrected_trimmed.tre")

#reduce to only polemoniaceae species
#non_phlox <- c(grep("Saurauia",phy$tip.label, value = TRUE),
#               grep("Actinidia",phy$tip.label, value = TRUE),
#               grep("Fouquieria",phy$tip.label, value = TRUE))
#phy <- drop.tip(phy, non_phlox)
#phy

#set colour palette
colours <- colorRampPalette(c("steelblue2", "paleturquoise3", 
                              "palegreen2", "yellow2", "salmon1", "darkorange", "red", 
                              "red4"))(100)

#get tip rates
tip_rates <- CladsOutput$lambdatip_map

#density plot of tip rates
car::densityPlot(tip_rates)

#get branch rates
rates <- CladsOutput$lambdai_map

#get colours for branch rates
col <- round((rates - min(rates))/diff(range(rates)) * 99) + 1

#plot phylogenetic tree with branches coloured by rates (should be same as Julia function 'plot_CladsOutput()')
png(file="figures/clads.png",width=15,height=15, units="in", res=200)
plot(phy, edge.color = colours[col], cex = 0.5, type="fan",edge.width=3,no.margin = TRUE)
dev.off()


#write output table
tip_rate_table <- data.frame(phy$tip.label,tip_rates)
colnames(tip_rate_table) <- c("names","lambda")
write.csv(tip_rate_table, "outputs/clads_tip_rate_table.csv")


