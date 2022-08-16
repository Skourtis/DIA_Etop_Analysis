#### CellROX Green Analysis ####
library(tidyverse)
library(data.table)
library(gridExtra)
# install.packages("interactions")
library(interactions)
list_of_significance <- list()
X0Hrs <- read_tsv(here::here("Datasets","Raw","220529_1 Meas5 Ev10 - shPRDX1 0hR.txt"), skip = 9) %>%
    dplyr::select(matches("uclei|Row|Column")) %>% 
    mutate(position = paste(Row,Column, sep = "_")) %>% janitor::clean_names()%>% 
    mutate(position = paste(row,column, sep = "_"),
           column = as.factor(column),
           Integrated_nuclear_488 = nuclei_selected_selected_nucleus_area_mm2*nuclei_selected_selected_intensity_nucleus_alexa_488_mean,
           Integrated_cytoplasmic_488 = nuclei_selected_selected_cytoplasm_area_mm2*nuclei_selected_selected_intensity_cytoplasm_alexa_488_mean,
           Integrated_cell_488 =nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_488_mean,
           Nuclear_vs_cyto_488 = Integrated_nuclear_488/Integrated_cytoplasmic_488,
           Integrated_nuclear_633 =nuclei_selected_selected_nucleus_area_mm2*nuclei_selected_selected_intensity_nucleus_alexa_633_mean ,
           Integrated_cytoplasmic_633 = nuclei_selected_selected_cytoplasm_area_mm2*nuclei_selected_selected_intensity_cytoplasm_alexa_633_mean,
           Integrated_cell_633 = nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_633_mean,
           Nuclear_vs_cyto_633 = Integrated_nuclear_633/Integrated_cytoplasmic_633
    ) 
X24Hrs <- read_tsv(here::here("Datasets","Raw","220529_2 Meas4 Ev6 - shPRDX1 24hR.txt"), skip = 9)%>%
    dplyr::select(matches("uclei|Row|Column"))%>% 
    mutate(position = paste(Row,Column, sep = "_"))%>% janitor::clean_names()%>% 
    mutate(position = paste(row,column, sep = "_"),
           column = as.factor(column),
           Integrated_nuclear_488 = nuclei_selected_selected_nucleus_area_mm2*nuclei_selected_selected_intensity_nucleus_alexa_488_mean,
           Integrated_cytoplasmic_488 = nuclei_selected_selected_cytoplasm_area_mm2*nuclei_selected_selected_intensity_cytoplasm_alexa_488_mean,
           Integrated_cell_488 =nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_488_mean,
           Nuclear_vs_cyto_488 = Integrated_nuclear_488/Integrated_cytoplasmic_488,
           Integrated_nuclear_633 =nuclei_selected_selected_nucleus_area_mm2*nuclei_selected_selected_intensity_nucleus_alexa_633_mean ,
           Integrated_cytoplasmic_633 = nuclei_selected_selected_cytoplasm_area_mm2*nuclei_selected_selected_intensity_cytoplasm_alexa_633_mean,
           Integrated_cell_633 = nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_633_mean,
           Nuclear_vs_cyto_633 = Integrated_nuclear_633/Integrated_cytoplasmic_633
    ) 
position <- openxlsx::read.xlsx(here::here("Datasets","Raw","220529 well labels.xlsx"))

#hist###
# samples <- c("220529_1 Meas5 - 0hR_2 - Ev6.txt","220529_2 Meas1 - 24hR_1 - Ev2.txt","220529_2 Meas4 - 24hR_2 - Ev3.txt","220529_1 Meas2 - 0hR_1 - Ev2.txt")
# for(i in samples){
# hist_df <- read_tsv(here::here("Datasets","Raw",i), skip = 9) %>%
#     dplyr::select(matches("uclei|Row|Column")) %>% 
#     mutate(position = paste(Row,Column, sep = "_")) %>% janitor::clean_names()
# hist_df$nuclei_selected_selected_intensity_nucleus_alexa_633_mean %>% hist(xlim=c(0, 50000), 
#                                                                            main = glue::glue("nucleus 633 {i}"), breaks = 100)
# hist_df$nuclei_selected_selected_intensity_cytoplasm_alexa_633_mean %>% hist(xlim=c(0, 50000), 
#                                                                              main = glue::glue("cyto 633 {i}"), breaks = 100)
# hist_df$nuclei_selected_selected_intensity_cell_alexa_633_mean %>% hist(xlim=c(0, 50000), 
#                                                                         main = glue::glue("cell 633 {i} "), breaks = 100)
# }



Combns <- data.table(X1 = "Comb1",
                     X2 = "Green_488",
                     X3 = "Mitotracker_633")
conditions <- data.frame(condition =  position[5:7,5:ncol(position)] %>% unlist(),
                         position =  position[8:10,5:ncol(position)] %>% unlist()) %>% 
    subset(str_detect(condition,"Comb1")) %>% 
    mutate(position = str_replace(position,"^4_","2_"))
X0Hrs <- left_join(X0Hrs,conditions) %>% na.omit() %>% separate(condition, sep = "_",into = c("Rep","sh","Treatment","Combo")) %>% 
    mutate(Treatment_bin = if_else(Treatment == 'DMSO',0,1),
           sh_bin = if_else(sh == "shNTC",0,1),
           sh_Treatment = paste(Treatment,sh, sep = "_") )
Combo_names_488 <- as_labeller(c(
    Comb1="Green_488"
))
Combo_names_633 <- as_labeller(c(
    Comb1="Mitotracker_633"
))

X0Hrs %>% 
    group_by(Combo,sh_Treatment) %>% 
    count()
ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_formula_nuc_cyt_488)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo", labeller = Combo_names_488, scales = "free_y")+ggtitle("T0  CellRox Green 488 Ratio N/C")
ggsave(here::here("Output","T0  CellRox Green 488.pdf"), width = 15, height = 11)
ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment,  y =  log2(nuclei_selected_selected_formula_nuc_cyt_633)))+
    
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo",labeller = Combo_names_633, scales= "free_y")+ggtitle("T0 633 Ratio N/C")
ggsave(here::here("Output","T0 633.pdf"), width = 15, height = 11)

#### nuclear and cyto seperately ####

ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(Integrated_nuclear_488)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo", labeller = Combo_names_488, scales = "free_y")+ggtitle("T0  CellRox Green 488 Only N")
ggsave(here::here("Output","T0  CellRox Green 488 Nuclear.pdf"), width = 15, height = 11)
ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(Integrated_cytoplasmic_488)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo", labeller = Combo_names_488, scales = "free_y")+ggtitle("T0  CellRox Green 488 Only C")
ggsave(here::here("Output","T0  CellRox Green 488 cyto.pdf"), width = 15, height = 11)


ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment,  y =  log2(Integrated_nuclear_633)))+
    
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo",labeller = Combo_names_633, scales= "free_y")+ggtitle("T0 633 Ratio  Only N")
ggsave(here::here("Output","T0 633  Only N.pdf"), width = 15, height = 11)

ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cytoplasm_alexa_488_mean)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo", labeller = Combo_names_488, scales = "free_y")+ggtitle("T0  CellRox Green 488 Only C")
ggsave(here::here("Output","T0  CellRox Green 488 Cyto.pdf"), width = 15, height = 11)
ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment,  y =  log2(nuclei_selected_selected_intensity_cytoplasm_alexa_633_mean)))+
    
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo",labeller = Combo_names_633, scales= "free_y")+ggtitle("T0 633 Ratio  Only C")
ggsave(here::here("Output","T0 633  Only C.pdf"), width = 15, height = 11)

#####


X24Hrs <- left_join(X24Hrs,conditions) %>% na.omit() %>% separate(condition, sep = "_",into = c("Rep","sh","Treatment","Combo")) %>% 
    mutate(Treatment_bin = if_else(Treatment == 'DMSO',0,1),
           sh_bin = if_else(sh == "shNTC",0,1),
           sh_Treatment = paste(Treatment,sh, sep = "_")) %>% na.omit()
ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(Integrated_nuclear_488)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    
    facet_wrap("Combo",scales = "free_y", labeller = Combo_names_488)+ggtitle("T24  CellRox Green 488 Ratio N/C")
ggsave(here::here("Output","T24 CellRox Green 488.pdf"), width = 15, height = 11)
ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_formula_nuc_cyt_633)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    
    facet_wrap("Combo",scales = "free_y", labeller = Combo_names_633)+ggtitle("T24 633 Ratio N/C")
ggsave(here::here("Output","T24 633.pdf"), width = 15, height = 11)

#### nuclear and cyto seperately ####

ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_nucleus_alexa_488_mean)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo", labeller = Combo_names_488, scales = "free_y")+ggtitle("T24  CellRox Green 488 Only N")
ggsave(here::here("Output","T24  CellRox Green 488 Nuclear.pdf"), width = 15, height = 11)
ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment,  y =  log2(nuclei_selected_selected_intensity_nucleus_alexa_633_mean)))+
    
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo",labeller = Combo_names_633, scales= "free_y")+ggtitle("T24 633 Ratio  Only N")
ggsave(here::here("Output","T24 633  Only N.pdf"), width = 15, height = 11)

ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cytoplasm_alexa_488_mean)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo", labeller = Combo_names_488, scales = "free_y")+ggtitle("T24  CellRox Green 488 Only C")
ggsave(here::here("Output","T24  CellRox Green 488 Cyto.pdf"), width = 15, height = 11)
ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment,  y =  log2(nuclei_selected_selected_intensity_cytoplasm_alexa_633_mean)))+
    
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo",labeller = Combo_names_633, scales= "free_y")+ggtitle("T24 633 Ratio  Only C")
ggsave(here::here("Output","T24 633  Only C.pdf"), width = 15, height = 11)

#####

ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cell_alexa_488_mean)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo", labeller = Combo_names_488, scales = "free_y")+ggtitle("T0 Green 488 Cell")
ggsave(here::here("Output","T0 Green Cell 488.pdf"), width = 15, height = 11)
ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment,  y =  log2(nuclei_selected_selected_intensity_cell_alexa_633_mean)))+
    
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo",labeller = Combo_names_633, scales= "free_y")+ggtitle("T0 633 Cell")
ggsave(here::here("Output","T0 Cell 633.pdf"), width = 15, height = 11)

ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cell_alexa_488_mean)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    
    facet_wrap("Combo",scales = "free_y", labeller = Combo_names_488)+ggtitle("T24 Green 488 Cell")
ggsave(here::here("Output","T24 Green cell 488.pdf"), width = 15, height = 11)

ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cell_alexa_633_mean)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    
    facet_wrap("Combo",scales = "free_y", labeller = Combo_names_633)+ggtitle("T24 633 Cell")
ggsave(here::here("Output","T24 Cell 633.pdf"), width = 15, height = 11)

Green_T0_T24 <- rbind(X0Hrs%>% subset(Combo == "Comb1")  %>% mutate(Time = "T0hrs"),
                        X24Hrs%>% subset(Combo == "Comb1")  %>% mutate(Time = "T24hrs"))
ggplot(Green_T0_T24 ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(Integrated_nuclear_488)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    
    facet_wrap("Time")+ggtitle("Green T0 T24 488 nuclear int")

ggsave(here::here("Output","Green T0 T24 nuclear int.pdf"), width = 7, height = 9)

ggplot(Green_T0_T24 ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(Integrated_cytoplasmic_488)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    
    facet_wrap("Time")+ggtitle("Green T0 T24 488 cyto int")

ggsave(here::here("Output","Green T0 T24 cyto int.pdf"), width = 7, height = 9)


ggplot(Green_T0_T24 ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(Integrated_cell_633)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    
    facet_wrap("Time")+ggtitle("Green T0 T24 633 nuclear int")

ggsave(here::here("Output","Green T0 T24 cell 633 int.pdf"), width = 7, height = 9)


ggplot(Green_T0_T24_nc ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(Integrated_nuclear_488)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    
    facet_wrap("Time")+ggtitle("Mtio  T0 T24 488 cell int")

ggsave(here::here("Output","Green T0 T24 nuclear int.pdf"), width = 7, height = 9)

tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     base_size = 10,
                     padding = unit(c(2, 4), "mm"))

for(i in unique(X0Hrs$Combo)){
    i  = "Comb1"
    Combination  = Combns[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
    model_0hrs_488 <- lm(log2(Integrated_nuclear_488)~sh_bin + Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X0Hrs %>% subset(
        Combo == i
    ) )
    # 
    interaction_plot_ <- interact_plot(model_0hrs_488, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination X0Hrs {Combination[1]}"))
    summary_tbl <- summary(model_0hrs_488)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X0Hrs Green {i} {Combination[1]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
    
    pdf(here::here("Output",glue::glue("Regression Green Images Comb {i}  X0Hrs {Combination[1]}.pdf")))
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
    model_0hrs_633 <- lm( nuclei_selected_selected_formula_nuc_cyt_633 ~sh_bin +Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X0Hrs %>% subset(
        Combo == i ))
    interaction_plot_ <- interact_plot(model_0hrs_633, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination  X0Hrs {Combination[2]}"))
    
    summary_tbl <- summary(model_0hrs_633)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X0Hrs {i} {Combination[2]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
    pdf(here::here("Output",glue::glue("Regression Images Comb X0Hrs {i} {Combination[2]}.pdf")))
    
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
}

for(i in unique(X24Hrs$Combo)){
    i  = "Comb1"
    Combination  = Combns[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
    model_24hrs_488 <- lm(log2(Integrated_nuclear_488)~sh_bin + Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X24Hrs %>% subset(
        Combo == i
    ) )
    # 
    interaction_plot_ <- interact_plot(model_24hrs_488, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination X24Hrs {Combination[1]}"))
    summary_tbl <- summary(model_24hrs_488)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X24Hrs Green {i} {Combination[1]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
    
    pdf(here::here("Output",glue::glue("Regression Green Images Comb {i}  X24Hrs {Combination[1]}.pdf")))
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
    model_0hrs_633 <- lm( nuclei_selected_selected_formula_nuc_cyt_633 ~sh_bin +Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X24Hrs %>% subset(
        Combo == i ))
    interaction_plot_ <- interact_plot(model_0hrs_633, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination  X24Hrs {Combination[2]}"))
    
    summary_tbl <- summary(model_0hrs_633)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X24Hrs {i} {Combination[2]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
    pdf(here::here("Output",glue::glue("Regression Images Comb X24Hrs {i} {Combination[2]}.pdf")))
    
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
}


for(i in unique(X0Hrs$Combo)){
    Combination  = Combns[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
    model_0hrs_488 <- lm(log2(nuclei_selected_selected_intensity_cell_alexa_488_mean)~sh_bin + Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X0Hrs %>% subset(
        Combo == i
    ) )
    # 
    interaction_plot_ <- interact_plot(model_0hrs_488, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination Cell X0Hrs {Combination[1]}"))
    summary_tbl <- summary(model_0hrs_488)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X0Hrs Green Cell {i} {Combination[1]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
    
    pdf(here::here("Output",glue::glue("Regression Green Cell Images Comb {i}  X0Hrs {Combination[1]}.pdf")))
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
    model_0hrs_633 <- lm( log2(nuclei_selected_selected_intensity_cell_alexa_633_mean) ~sh_bin +Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X0Hrs %>% subset(
        Combo == i ))
    interaction_plot_ <- interact_plot(model_0hrs_633, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination Cell X0Hrs {Combination[2]}"))
    
    summary_tbl <- summary(model_0hrs_633)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X0Hrs Cell {i} {Combination[2]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
    pdf(here::here("Output",glue::glue("Regression Images Comb Cell X0Hrs {i} {Combination[2]}.pdf")))
    
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
}

for(i in unique(X24Hrs$Combo)){
    Combination  = Combns[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
    model_24hrs_488 <- lm(log2(nuclei_selected_selected_intensity_cell_alexa_488_mean)~sh_bin + Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X24Hrs %>% subset(
        Combo == i
    ) )
    # 
    interaction_plot_ <- interact_plot(model_24hrs_488, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination X24Hrs Cell {Combination[1]}"))
    summary_tbl <- summary(model_24hrs_488)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X24Hrs Green Cell {i} {Combination[1]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
    
    pdf(here::here("Output",glue::glue("Regression Green Cell Images Comb {i}  X24Hrs {Combination[1]}.pdf")))
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
    model_0hrs_633 <- lm( log2(nuclei_selected_selected_intensity_cell_alexa_633_mean) ~sh_bin +Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X24Hrs %>% subset(
        Combo == i ))
    interaction_plot_ <- interact_plot(model_0hrs_633, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination  X24Hrs Cell {Combination[2]}"))
    
    summary_tbl <- summary(model_0hrs_633)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X24Hrs Cell {i} {Combination[2]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
    pdf(here::here("Output",glue::glue("Regression Images Comb X24Hrs Cell {i} {Combination[2]}.pdf")))
    
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
}

