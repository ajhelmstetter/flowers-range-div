# load libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(tidyverse)
library(ggpubr)
library(rstatix)

# Copy cleaned_trees.rds from tree_and_traits repo
# Copy species_with_ranges.rds from gbif-bulk repo

# Read in trees
cleaned_trees <- readRDS("data/raw-data/cleaned_trees.rds")

# Read in ranges
ranges <- readRDS("data/raw-data/species_with_ranges.rds")

# Read in trait data
#ss1 <- read.csv("../trees_and_traits/data/raw-data/traits/wang_et_al_2021_sexual_systems/SplistWithSexSyst.csv")
ss <- readRDS("../trees_and_traits/data/derived-data/cleaned_traits/sexual_system-Wang_et_al-2021.rds")

# Get overlapping species, trim, sort and combine
ranges <- ranges[ranges$species_name%in%ss$species_name,]
ss <- ss[ss$species_name%in%ranges$species_name,]

ss <- ss[order(ss$species_name),]
ranges <- ranges[order(ranges$species_name),]

table(ss$species_name==ranges$species_name)

ranges_ss <- cbind(ranges,ss)

# log range size
ranges_ss$log_range_size <- log(ranges_ss$range_size)

ranges_ss <- ranges_ss %>%
  reorder_levels(SexSyst, order = c("hermaphroditism", "monoecy", "dioecy"))

####
# ---- Violin plot ----
####

# https://r-graph-gallery.com/violin_and_boxplot_ggplot2.html

# sample size
sample_size = ranges_ss %>% group_by(SexSyst) %>% summarize(num=n())

# Plot
ranges_ss %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(SexSyst, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=log_range_size, fill=SexSyst)) +
  geom_violin(width=1) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")

####
# ---- Testing assumptions ----
####

# FROM: https://www.datanovia.com/en/lessons/anova-in-r/

# Build the linear model
model  <- lm(log_range_size ~ SexSyst, data = ranges_ss)

# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
# NOTE: Limited to 5000 - randomly sample 2500 points?
ranges_ss_2500 <- ranges_ss[sample(nrow(ranges_ss), 2500), ]
model_2500  <- lm(log_range_size ~ SexSyst, data = ranges_ss_2500)
ggqqplot(residuals(model_2500))
shapiro_test(residuals(model_2500))

# The residuals versus fits plot can be used to check the homogeneity of variances.
plot(model, 1)
# Levene’s test to check the homogeneity of variances:
ranges_ss %>% levene_test(log_range_size ~ SexSyst)

####
# ---- Kruskall-Wallis ----
####

# Non-parametric ANOVA
res.kw <- ranges_ss %>% kruskal_test(log_range_size ~ SexSyst)
res.kw

# Non-parametric tukey aka Dunn test 
pwc <- ranges_ss %>% dunn_test(log_range_size ~ SexSyst)
pwc

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "SexSyst")

ggboxplot(ranges_ss, x = "SexSyst", y = "log_range_size") +
  stat_pvalue_manual(pwc, hide.ns = FALSE, step.increase = 0.075, tip.length = 0.02) +
  labs(
    subtitle = get_test_label(res.kw, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

# Could improve plot with : https://r-graph-gallery.com/web-violinplot-with-ggstatsplot.html

####
# ---- ANOVA ----
####

# NOTE: Currently variances are not homogeneous so not run
# res.aov <- ranges_ss %>% anova_test(log_range_size ~ SexSyst)
# res.aov

# NOTE: Currently data is not normally distributed so not run
# # In a situation where the homogeneity of variance assumption is not met, 
# # you can compute the Welch one-way ANOVA test using the function welch_anova_test()
# res.aov <- ranges_ss %>% welch_anova_test(log_range_size ~ SexSyst)
# res.aov

####
# ---- Pairwise comparisons ----
####

# NOTE: Currently data is not normally distributed so not run
# # A significant one-way ANOVA is generally followed up by Tukey post-hoc tests 
# # to perform multiple pairwise comparisons between groups.
# pwc <- ranges_ss %>% tukey_hsd(log_range_size ~ SexSyst)
# pwc

# The output contains the following columns:
# estimate: estimate of the difference between means of the two groups
# conf.low, conf.high: the lower and the upper end point of the confidence interval at 95% (default)
# p.adj: p-value after adjustment for the multiple comparisons.