### QC screen #### 
if(!"MAGeCKFlute" %in% installed.packages()) BiocManager::install("MAGeCKFlute")
library(MAGeCKFlute)
library(tidyverse)
ouput_folder <- here::here("Output","Figures")
countsummary <- read.delim(here::here("Datasets","Processed","Amandine_screen","all.countsummary_screen1.txt"), check.names = FALSE)
countsummary <- countsummary %>% dplyr::mutate(Label = if_else(str_detect(File,"top5"),"H2gAX high","Unsorted"))
countsummary$Unmapped = countsummary$Reads - countsummary$Mapped
# gg = reshape2::melt(countsummary[, c("Label", "Mapped", "Unmapped")], id.vars = "Label")
# gg$variable = factor(gg$variable, levels = c("Unmapped", "Mapped"))
# gg = gg[order(gg$Label, gg$variable), ]
# p = BarView(gg, x = "Label", y = "value", fill = "variable", 
#             position = "stack", xlab = NULL, ylab = "Reads", main = "Map ratio")
# p + scale_fill_manual(values = c("#9BC7E9", "#1C6DAB"))
MapRatesView(countsummary)+ggtitle("Mapping Ratio H2gAX - Etoposide Treatment")
ggsave(here::here(ouput_folder,"Mapping Ratio H2gAX.pdf"))

countsummary$Missed = log10(countsummary$Zerocounts)
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")
ggsave(here::here(ouput_folder,"Missingness guides H2gAX.pdf"))

countsummary <- read.delim(here::here("Datasets","Processed","Amandine_screen","all.countsummary_screen2.txt"), check.names = FALSE)
countsummary <- countsummary %>% dplyr::mutate(Label = case_when(
    str_detect(File,"library")~ "Library",
    str_detect(File,"etop")~ "Etop 5d Recovery",
    str_detect(File,"S52965")~ "DMSO 5d Control")) %>% 
    mutate(Label = factor(Label, levels =  c("Library","DMSO 5d Control","Etop 5d Recovery")))
countsummary$Unmapped = countsummary$Reads - countsummary$Mapped
# gg = reshape2::melt(countsummary[, c("Label", "Mapped", "Unmapped")], id.vars = "Label")
# gg$variable = factor(gg$variable, levels = c("Unmapped", "Mapped"))
# gg = gg[order(gg$Label, gg$variable), ]
# p = BarView(gg, x = "Label", y = "value", fill = "variable", 
#             position = "stack", xlab = NULL, ylab = "Reads", main = "Map ratio")
# p + scale_fill_manual(values = c("#9BC7E9", "#1C6DAB"))
MapRatesView(countsummary)+ggtitle("Mapping Ratio Etop Release")
    
ggsave(here::here(ouput_folder,"Mapping Ratio D14.pdf"))

countsummary$Missed = log10(countsummary$Zerocounts)
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")
ggsave(here::here(ouput_folder,"Missingness guides D14.pdf"))



counts <- read.delim(here::here("Datasets","Processed","Amandine_screen","Amandine.count.txt"), check.names = FALSE)
counts %>% ggplot(aes(x= log10(`10d-etop`)))+
    geom_histogram(colour = "grey20",
                   fill = "grey50")+
    theme_bw(base_size = 25)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    labs(x = "Log10 Counts")+
    ggtitle("Unsorted Counts")
ggsave(here::here(ouput_folder,"Counts_hist_unsort.pdf"))

counts %>% ggplot(aes(x= log10(library)))+
    geom_histogram(colour = "grey20",
                   fill = "grey50")+
    
    theme_bw(base_size = 25)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    labs(x = "Log10 Counts")+
    ggtitle("Library Counts")
ggsave(here::here(ouput_folder,"Counts_hist_library.pdf"))
