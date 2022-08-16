library(tidyverse)
library(rstatix)
set.seed(12345)
#create artificial datasets of cells with number of foci
# 3 replicates for a cntrl and treatment conditions

Foci_per_cell_df <- data.frame(Condition = rep(c("Cntrl","Treatment"), each = 300),
           Replicate = rep(c(1,2,3), times = 200),
           Foci_per_cell = c(rpois(300,1), # control conditions are modelled as poisson distribution with lambda 1
                             rpois(300,3)))# while Treatment with lambda 2
Foci_per_cell_df %>% subset(Condition  == 'Cntrl') %>% pull(Foci_per_cell) %>% hist(main = "Control distribution")
Foci_per_cell_df %>% subset(Condition  == 'Treatment') %>% pull(Foci_per_cell) %>% hist(main = 'Treatment Distribution')
top5_perc_control <- Foci_per_cell_df %>% 
    subset(Condition == "Cntrl") %>% pull(Foci_per_cell) %>% 
    quantile(prob = seq(0, 1, length = 21), type = 5) %>% .[20] # the top 5% is in the 20th position
#all values above 3 (in this example) are considered positive
Foci_per_cell_df <- Foci_per_cell_df %>% mutate(State = if_else(Foci_per_cell>top5_perc_control, 1,0)) # 1 = positive, 0 = Negative
Foci_per_cell_df_sum <- Foci_per_cell_df %>% group_by(Replicate,Condition) %>% summarise(Total_positives = mean(State))
Foci_per_cell_df_sum %>% ungroup %>% 
    group_by(Condition) %>%
    shapiro_test(Total_positives) # pval needs to be above 0.05 to allow for test
Foci_per_cell_df_sum %>% ungroup() %>% levene_test(Total_positives ~ Condition) #needs to be more than 0.05

#after checking assumptions variance and normality we can perform a parametric t test
#if not we do a non-parametric test
t.test(Foci_per_cell_df_sum %>% subset(Condition  == 'Cntrl') %>% pull(Total_positives), 
       Foci_per_cell_df_sum %>% subset(Condition  == 'Treatment') %>% pull(Total_positives))
