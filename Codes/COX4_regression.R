#### COX4 Analysis ####
library(tidyverse)
library(data.table)
library(gridExtra)
# install.packages("interactions")
library(interactions)
list_of_significance <- list()
X0Hrs <- read_tsv(here::here("Datasets","Raw","220529_1 0hR Meas1 Ev6 - confocal.txt"), skip = 9) %>%
    dplyr::select(matches("uclei|Row|Column")) %>% 
    mutate(position = paste(Row,Column, sep = "_")) %>% janitor::clean_names()
X0Hrs <- X0Hrs %>% 
  mutate(nuclei_selected_selected_intensity_nucleus_alexa_488_int = nuclei_selected_selected_intensity_nucleus_alexa_488_mean *nuclei_selected_selected_nucleus_area_mm2,
         nuclei_selected_selected_intensity_cytoplasm_alexa_488_int  = nuclei_selected_selected_intensity_cytoplasm_alexa_488_mean*nuclei_selected_selected_cytoplasm_area_mm2,
         nuclei_selected_selected_intensity_cell_alexa_488_int= nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_488_mean,
          nuclei_selected_selected_formula_nuc_cyt_488 =nuclei_selected_selected_intensity_nucleus_alexa_488_int/ nuclei_selected_selected_intensity_cytoplasm_alexa_488_int,
         nuclei_selected_selected_intensity_nucleus_alexa_568_int = nuclei_selected_selected_intensity_nucleus_alexa_568_mean *nuclei_selected_selected_nucleus_area_mm2,
         nuclei_selected_selected_intensity_cytoplasm_alexa_568_int  = nuclei_selected_selected_intensity_cytoplasm_alexa_568_mean*nuclei_selected_selected_cytoplasm_area_mm2,
         nuclei_selected_selected_intensity_cell_alexa_568_int= nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_568_mean,
         nuclei_selected_selected_formula_nuc_cyt_568 =nuclei_selected_selected_intensity_nucleus_alexa_568_int/ nuclei_selected_selected_intensity_cytoplasm_alexa_568_int
  )
X24Hrs <- read_tsv(here::here("Datasets","Raw","220529_2 24hR Meas1 Ev4 - confocal.txt"), skip = 9)%>%
    dplyr::select(matches("uclei|Row|Column"))%>% 
    mutate(position = paste(Row,Column, sep = "_"))%>% janitor::clean_names()%>% 
  mutate(nuclei_selected_selected_intensity_nucleus_alexa_488_int = nuclei_selected_selected_intensity_nucleus_alexa_488_mean *nuclei_selected_selected_nucleus_area_mm2,
         nuclei_selected_selected_intensity_cytoplasm_alexa_488_int  = nuclei_selected_selected_intensity_cytoplasm_alexa_488_mean*nuclei_selected_selected_cytoplasm_area_mm2,
         nuclei_selected_selected_intensity_cell_alexa_488_int= nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_488_mean,
         nuclei_selected_selected_formula_nuc_cyt_488 =nuclei_selected_selected_intensity_nucleus_alexa_488_int/ nuclei_selected_selected_intensity_cytoplasm_alexa_488_int,
         nuclei_selected_selected_intensity_nucleus_alexa_568_int = nuclei_selected_selected_intensity_nucleus_alexa_568_mean *nuclei_selected_selected_nucleus_area_mm2,
         nuclei_selected_selected_intensity_cytoplasm_alexa_568_int  = nuclei_selected_selected_intensity_cytoplasm_alexa_568_mean*nuclei_selected_selected_cytoplasm_area_mm2,
         nuclei_selected_selected_intensity_cell_alexa_568_int= nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_568_mean,
         nuclei_selected_selected_formula_nuc_cyt_568 =nuclei_selected_selected_intensity_nucleus_alexa_568_int/ nuclei_selected_selected_intensity_cytoplasm_alexa_568_int
  )

position <- openxlsx::read.xlsx(here::here("Datasets","Raw","220529 well labels.xlsx"))

Combns <- position[5:7,1:3] %>% as.data.table()
conditions <- data.frame(condition =  position[5:7,5:ncol(position)] %>% unlist(),
           position =  position[8:10,5:ncol(position)] %>% unlist())



X0Hrs <- left_join(X0Hrs,conditions) %>% na.omit() %>% separate(condition, sep = "_",into = c("Rep","sh","Treatment","Combo")) %>% 
    mutate(Treatment_bin = if_else(Treatment == 'DMSO',0,1),
           sh_bin = if_else(sh == "shNTC",0,1),
           sh_Treatment = paste(Treatment,sh, sep = "_") )
X24Hrs <- left_join(X24Hrs,conditions) %>% separate(condition, sep = "_",into = c("Rep","sh","Treatment","Combo")) %>% 
  mutate(Treatment_bin = if_else(Treatment == 'DMSO',0,1),
         sh_bin = if_else(sh == "shNTC",0,1),
         sh_Treatment = paste(Treatment,sh, sep = "_")) %>% na.omit()
gh2ax_signal <- X0Hrs %>% subset(Combo == "Comb1") %>% pull(nuclei_selected_selected_intensity_nucleus_alexa_568_int) %>% quantile(type = 5, probs = seq(0,1,0.1))
X0Hrs_gh2ax <- X0Hrs %>% 
  mutate(gh2ax_level = case_when(
    nuclei_selected_selected_intensity_nucleus_alexa_568_int<gh2ax_signal[4]~"low gh2ax",
    nuclei_selected_selected_intensity_nucleus_alexa_568_int>gh2ax_signal[8]~"high gh2ax",
    TRUE~"middle signal"))

gh2ax_signal <- X24Hrs %>% subset(Combo == "Comb1") %>% pull(nuclei_selected_selected_intensity_nucleus_alexa_568_int) %>% quantile(type = 5, probs = seq(0,1,0.1))
X24Hrs_gh2ax <- X24Hrs %>% 
  mutate(gh2ax_level = case_when(
    nuclei_selected_selected_intensity_nucleus_alexa_568_int<gh2ax_signal[4]~"low gh2ax",
    nuclei_selected_selected_intensity_nucleus_alexa_568_int>gh2ax_signal[8]~"high gh2ax",
    TRUE~"middle signal"))

#### testing#####
shNTC_T0 <- conditions %>% subset(str_detect(condition,"shNTC")& str_detect(condition,"Comb1") )
X0Hrs %>% subset(position %in% shNTC_T0$position) %>% 
  ggplot(aes(x = position, colour = Treatment,y = log2(nuclei_selected_selected_formula_nuc_cyt_488)))+
  geom_boxplot()+
  ggtitle("T0 N/C shCTRL")

X24Hrs %>% subset(position %in% shNTC_T0$position) %>% 
  ggplot(aes(x = position, colour = Treatment,y = log2(nuclei_selected_selected_formula_nuc_cyt_488)))+
  geom_boxplot()+
  ggtitle("T24 N/C shCTRL")
####

Combo_names_488 <- as_labeller(c(
    Comb1="PRDX1",
    Comb2="COX4",
    Comb3="COX4"
))
Combo_names_568 <- as_labeller(c(
    Comb1="gH2AX",
    Comb2="gH2AX",
    Comb3="SDHB"
))
ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_formula_nuc_cyt_488)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                 "shPRDX1" = "black"))+
    facet_wrap("Combo", labeller = Combo_names_488, scales = "free_y")+ggtitle("T0 488 Ratio N/C")
ggsave(here::here("Output","T0 488_int.pdf"), width = 15, height = 11)
# ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment,  y =  log2(nuclei_selected_selected_formula_nuc_cyt_568)))+
#   
#     geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
#                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#     scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
#                                  "Etoposide" = scales::muted("red")))+
#     scale_colour_manual(values = c("shNTC" ="grey20" ,
#                                    "shPRDX1" = "black"))+
#     facet_wrap("Combo",labeller = Combo_names_568, scales= "free_y")+ggtitle("T0 568 Ratio N/C")
# ggsave(here::here("Output","T0 568_int.pdf"), width = 15, height = 11)
    
ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_formula_nuc_cyt_488)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    
    facet_wrap("Combo",scales = "free_y", labeller = Combo_names_488)+ggtitle("T24 488 Ratio N/C")
ggsave(here::here("Output","T24 488_int.pdf"), width = 15, height = 11)
ggplot(X24Hrs_gh2ax %>% subset(Combo == "Comb1" & (gh2ax_level !="middle signal") & str_detect(sh_Treatment,"NTC")) ,
       aes(x =gh2ax_level, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_formula_nuc_cyt_488)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  scale_colour_manual(values = c("shNTC" ="grey20" ,
                                 "shPRDX1" = "black"))+
  
  facet_wrap("sh_Treatment")+ggtitle("T24 488 Ratio N/C")
ggsave(here::here("Output","T24 nc 488_int highgh2ax.pdf"), width = 4, height = 6)

ggplot(X24Hrs_gh2ax %>% subset(Combo == "Comb1") ,
       aes(x =gh2ax_level, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_nucleus_alexa_568_int)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  scale_colour_manual(values = c("shNTC" ="grey20" ,
                                 "shPRDX1" = "black"))+
  
  facet_wrap("sh_Treatment")+ggtitle("T24 568 nuclear_int")

ggplot(X0Hrs_gh2ax %>% subset(Combo == "Comb1" & (gh2ax_level !="middle signal") & str_detect(sh_Treatment,"NTC")) ,
       aes(x =gh2ax_level, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_formula_nuc_cyt_488)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  scale_colour_manual(values = c("shNTC" ="grey20" ,
                                 "shPRDX1" = "black"))+
  
  facet_wrap("sh_Treatment")+ggtitle("T0 488 Ratio N/C")
ggsave(here::here("Output","T0 nc 488_int highgh2ax.pdf"), width = 4, height = 6)

ggplot(X0Hrs_gh2ax %>% subset(Combo == "Comb1") ,
       aes(x =gh2ax_level, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_nucleus_alexa_568_int)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  scale_colour_manual(values = c("shNTC" ="grey20" ,
                                 "shPRDX1" = "black"))+
  
  facet_wrap("sh_Treatment")+ggtitle("T0 568_nuclear_int")


# ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_formula_nuc_cyt_568)))+
#     geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
#                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#     scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
#                                  "Etoposide" = scales::muted("red")))+
#     scale_colour_manual(values = c("shNTC" ="grey20" ,
#                                    "shPRDX1" = "black"))+
#     
#     facet_wrap("Combo",scales = "free_y", labeller = Combo_names_568)+ggtitle("T24 568 Ratio N/C")
# ggsave(here::here("Output","T24 568_int.pdf"), width = 15, height = 11)

ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cell_alexa_488_int)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo", labeller = Combo_names_488, scales = "free_y")+ggtitle("T0 488 Cell")
ggsave(here::here("Output","T0 Cell 488_int.pdf"), width = 15, height = 11)
ggplot(X0Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment,  y =  log2(nuclei_selected_selected_intensity_cell_alexa_568_int)))+
    
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    facet_wrap("Combo",labeller = Combo_names_568, scales= "free_y")+ggtitle("T0 568 Cell")
ggsave(here::here("Output","T0 Cell 568_int.pdf"), width = 15, height = 11)

ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cell_alexa_488_int)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    
    facet_wrap("Combo",scales = "free_y", labeller = Combo_names_488)+ggtitle("T24 488 Cell")
ggsave(here::here("Output","T24 cell 488_int.pdf"), width = 15, height = 11)
ggplot(X24Hrs ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cell_alexa_568_int)))+
    geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                                 "Etoposide" = scales::muted("red")))+
    scale_colour_manual(values = c("shNTC" ="grey20" ,
                                   "shPRDX1" = "black"))+
    
    facet_wrap("Combo",scales = "free_y", labeller = Combo_names_568)+ggtitle("T24 568 Cell")
ggsave(here::here("Output","T24 Cell 568_int.pdf"), width = 15, height = 11)

PRDX1_T0_T24_nc <- rbind(X0Hrs%>% subset(Combo == "Comb1")  %>% mutate(Time = "T0hrs"),
                        X24Hrs%>% subset(Combo == "Comb1")  %>% mutate(Time = "T24hrs"))
ggplot(PRDX1_T0_T24_nc ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cell_alexa_488_int)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  scale_colour_manual(values = c("shNTC" ="grey20" ,
                                 "shPRDX1" = "black"))+
  
  facet_wrap("Time")+ggtitle("PRDX1 T0 T24 488 Cell")
COX4_T0_T24_nc <- rbind(X0Hrs%>% subset(Combo == "Comb2")  %>% mutate(Time = "T0hrs"),
                        X24Hrs%>% subset(Combo == "Comb2")  %>% mutate(Time = "T24hrs"))

ggplot(COX4_T0_T24_nc ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_nucleus_alexa_488_int)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  scale_colour_manual(values = c("shNTC" ="grey20" ,
                                 "shPRDX1" = "black"))+
  
  facet_wrap("Time")+ggtitle("Cox4 T0 T24 488 Nuclear")

ggsave(here::here("Output","COX4 T0 T24 nuclear_int.pdf"), width = 7, height = 9)

ggplot(COX4_T0_T24_nc ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cytoplasm_alexa_488_int)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  scale_colour_manual(values = c("shNTC" ="grey20" ,
                                 "shPRDX1" = "black"))+
  
  facet_wrap("Time")+ggtitle("Cox4 T0 T24 488 cyto")

ggsave(here::here("Output","COX4 T0 T24 cyto_int.pdf"), width = 7, height = 9)


gH2ax_T0_T24_nc <- rbind(X0Hrs%>% subset(Combo == "Comb2")  %>% mutate(Time = "T0hrs"),
                        X24Hrs%>% subset(Combo == "Comb2")  %>% mutate(Time = "T24hrs"))
ggplot(gH2ax_T0_T24_nc ,aes(x =sh_Treatment, colour = sh, fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_nucleus_alexa_568_int)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  scale_colour_manual(values = c("shNTC" ="grey20" ,
                                 "shPRDX1" = "black"))+
  
  facet_wrap("Time")+ggtitle("gH2ax T0 T24 568 Nucleus")

ggsave(here::here("Output","gH2ax T0 T24 nuclear_int.pdf"), width = 7, height = 9)

tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     base_size = 10,
                     padding = unit(c(2, 4), "mm"))

for(i in unique(X0Hrs$Combo)){
  i = unique(X0Hrs$Combo)[2]
Combination  = Combns[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
model_0hrs_488 <- lm(log2(nuclei_selected_selected_intensity_nucleus_alexa_488_int)~sh_bin + Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X0Hrs %>% subset(
    Combo == i
) )
 # 
interaction_plot_ <- interact_plot(model_0hrs_488, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination X0Hrs {Combination[1]}"))
summary_tbl <- summary(model_0hrs_488)
summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
list_of_significance[[glue::glue("X0Hrs {i} {Combination[1]}")]] <- summary_tbl
tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)

pdf(here::here("Output",glue::glue("Regression Images Comb {i}  X0Hrs {Combination[1]}_int.pdf")))
grid.arrange(interaction_plot_, tbl, 
             nrow = 2, heights = c(2, 2))
dev.off()
model_0hrs_568 <- lm( nuclei_selected_selected_formula_nuc_cyt_568 ~sh_bin +Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X0Hrs %>% subset(
    Combo == i ))
interaction_plot_ <- interact_plot(model_0hrs_568, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination X0Hrs {Combination[2]}"))

summary_tbl <- summary(model_0hrs_568)
summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
list_of_significance[[glue::glue("X0Hrs {i} {Combination[2]}")]] <- summary_tbl
tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
pdf(here::here("Output",glue::glue("Regression Images Comb X0Hrs {i} {Combination[2]}_int.pdf")))

grid.arrange(interaction_plot_, tbl, 
             nrow = 2, heights = c(2, 2))
dev.off()
}

for(i in na.omit(unique(X24Hrs$Combo))){
    Combination  = Combns[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
    model_0hrs_488 <- lm(nuclei_selected_selected_formula_nuc_cyt_488~sh_bin + Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X24Hrs %>% subset(
        Combo == i
    ) )
    # 
    interaction_plot_ <- interact_plot(model_0hrs_488, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination X24Hrs {Combination[1]}"))
    summary_tbl <- summary(model_0hrs_488)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
    list_of_significance[[glue::glue("X24Hrs {i} {Combination[1]}")]] <- summary_tbl
   pdf(here::here("Output",glue::glue("Regression Images Comb {i} X24Hrs {Combination[1]}_int.pdf")))
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
    model_0hrs_568 <- lm( nuclei_selected_selected_formula_nuc_cyt_568 ~sh_bin +Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X24Hrs %>% subset(
        Combo == i ))
    interaction_plot_ <- interact_plot(model_0hrs_568, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination X24Hrs {Combination[2]}"))
    
    summary_tbl <- summary(model_0hrs_568)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X24Hrs {i} {Combination[2]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
   pdf(here::here("Output",glue::glue("Regression Images Comb X24Hrs {i} {Combination[2]}_int.pdf")))
    
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
}

for(i in unique(X0Hrs$Combo)){
    Combination  = Combns[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
    model_0hrs_488 <- lm(log2(nuclei_selected_selected_intensity_cell_alexa_488_mean)~sh_bin + Rep+ Treatment_bin + Treatment_bin*sh_bin,
                         data =X0Hrs %>% subset(
        Combo == i
    ) )
    # 
    interaction_plot_ <- interact_plot(model_0hrs_488, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination X0Hrs Cell {Combination[1]}"))
    summary_tbl <- summary(model_0hrs_488)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X0Hrs Cell {i} {Combination[1]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
    
   pdf(here::here("Output",glue::glue("Regression Images Comb {i}  X0Hrs Cell {Combination[1]}_int.pdf")))
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
    model_0hrs_568 <- lm(log2(nuclei_selected_selected_intensity_cell_alexa_568_mean) ~sh_bin +Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X0Hrs %>% subset(
        Combo == i ))
    interaction_plot_ <- interact_plot(model_0hrs_568, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination X0Hrs Cell {Combination[2]}"))
    
    summary_tbl <- summary(model_0hrs_568)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X0Hrs Cell {i} {Combination[2]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
   pdf(here::here("Output",glue::glue("Regression Images Comb X0Hrs {i} Cell {Combination[2]}_int.pdf")))
    
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
}

for(i in na.omit(unique(X24Hrs$Combo))){
    Combination  = Combns[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
    model_0hrs_488 <- lm(log2(nuclei_selected_selected_intensity_cell_alexa_488_int)~sh_bin + Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X24Hrs %>% subset(
        Combo == i
    ) )
    # 
    interaction_plot_ <- interact_plot(model_0hrs_488, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination X24Hrs Cell {Combination[1]}"))
    summary_tbl <- summary(model_0hrs_488)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
    list_of_significance[[glue::glue("X24Hrs Cell {i} {Combination[1]}")]] <- summary_tbl
   pdf(here::here("Output",glue::glue("Regression Images Comb {i} Cell X24Hrs {Combination[1]}_int.pdf")))
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
    model_0hrs_568 <- lm( log2(nuclei_selected_selected_intensity_cell_alexa_568_mean) ~sh_bin +Rep+ Treatment_bin + Treatment_bin*sh_bin, data =X24Hrs %>% subset(
        Combo == i ))
    interaction_plot_ <- interact_plot(model_0hrs_568, pred = sh_bin, modx = Treatment_bin)+ggtitle(glue::glue("Combination X24Hrs Cell {Combination[2]}"))
    
    summary_tbl <- summary(model_0hrs_568)
    summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
    list_of_significance[[glue::glue("X24Hrs Cell {i} {Combination[2]}")]] <- summary_tbl
    tbl <- tableGrob(summary_tbl, rows=NULL, theme=tt)
   pdf(here::here("Output",glue::glue("Regression Images Comb X24Hrs {i} Cell {Combination[2]}_int.pdf")))
    
    grid.arrange(interaction_plot_, tbl, 
                 nrow = 2, heights = c(2, 2))
    dev.off()
}

###### wt ##### 
X0Hrs_wt <- read_tsv(here::here("Datasets","Raw","220529_1 0hR Meas1 Ev6 - confocal.txt"), skip = 9) %>%
  dplyr::select(matches("uclei|Row|Column")) %>% 
  mutate(position_wt = paste(Row,Column, sep = "_")) %>% janitor::clean_names()
X0Hrs_wt <- X0Hrs_wt %>% 
  mutate(nuclei_selected_selected_intensity_nucleus_alexa_488_int = nuclei_selected_selected_intensity_nucleus_alexa_488_mean *nuclei_selected_selected_nucleus_area_mm2,
         nuclei_selected_selected_intensity_cytoplasm_alexa_488_int  = nuclei_selected_selected_intensity_cytoplasm_alexa_488_mean*nuclei_selected_selected_cytoplasm_area_mm2,
         nuclei_selected_selected_intensity_cell_alexa_488_int= nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_488_mean,
         nuclei_selected_selected_formula_nuc_cyt_488 =nuclei_selected_selected_intensity_nucleus_alexa_488_int/ nuclei_selected_selected_intensity_cytoplasm_alexa_488_int,
         nuclei_selected_selected_intensity_nucleus_alexa_568_int = nuclei_selected_selected_intensity_nucleus_alexa_568_mean *nuclei_selected_selected_nucleus_area_mm2,
         nuclei_selected_selected_intensity_cytoplasm_alexa_568_int  = nuclei_selected_selected_intensity_cytoplasm_alexa_568_mean*nuclei_selected_selected_cytoplasm_area_mm2,
         nuclei_selected_selected_intensity_cell_alexa_568_int= nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_568_mean,
         nuclei_selected_selected_formula_nuc_cyt_568 =nuclei_selected_selected_intensity_nucleus_alexa_568_int/ nuclei_selected_selected_intensity_cytoplasm_alexa_568_int
  )
X24Hrs_wt <- read_tsv(here::here("Datasets","Raw","220529_2 24hR Meas1 Ev4 - confocal.txt"), skip = 9)%>%
  dplyr::select(matches("uclei|Row|Column"))%>% 
  mutate(position_wt = paste(Row,Column, sep = "_"))%>% janitor::clean_names()%>% 
  mutate(nuclei_selected_selected_intensity_nucleus_alexa_488_int = nuclei_selected_selected_intensity_nucleus_alexa_488_mean *nuclei_selected_selected_nucleus_area_mm2,
         nuclei_selected_selected_intensity_cytoplasm_alexa_488_int  = nuclei_selected_selected_intensity_cytoplasm_alexa_488_mean*nuclei_selected_selected_cytoplasm_area_mm2,
         nuclei_selected_selected_intensity_cell_alexa_488_int= nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_488_mean,
         nuclei_selected_selected_formula_nuc_cyt_488 =nuclei_selected_selected_intensity_nucleus_alexa_488_int/ nuclei_selected_selected_intensity_cytoplasm_alexa_488_int,
         nuclei_selected_selected_intensity_nucleus_alexa_568_int = nuclei_selected_selected_intensity_nucleus_alexa_568_mean *nuclei_selected_selected_nucleus_area_mm2,
         nuclei_selected_selected_intensity_cytoplasm_alexa_568_int  = nuclei_selected_selected_intensity_cytoplasm_alexa_568_mean*nuclei_selected_selected_cytoplasm_area_mm2,
         nuclei_selected_selected_intensity_cell_alexa_568_int= nuclei_selected_selected_cell_area_mm2*nuclei_selected_selected_intensity_cell_alexa_568_mean,
         nuclei_selected_selected_formula_nuc_cyt_568 =nuclei_selected_selected_intensity_nucleus_alexa_568_int/ nuclei_selected_selected_intensity_cytoplasm_alexa_568_int
  )
position_wt <- openxlsx::read.xlsx(here::here("Datasets","Raw","220529 well labels.xlsx"))
Combns_wt <- position_wt[17:18,1:3] %>% as.data.table()
conditions <- data.frame(condition =  position_wt[17:18,5:ncol(position_wt)] %>% unlist(),
                         position_wt =  position_wt[19:20,5:ncol(position_wt)] %>% unlist())
X0Hrs_wt <- left_join(X0Hrs_wt,conditions) %>% na.omit() %>% separate(condition, sep = "_",into = c("Rep","Trep","WT","Treatment","Combo")) %>% 
  mutate(Treatment_bin = if_else(Treatment == 'DMSO',0,1)) 
Combo_names_488 <- as_labeller(c(
  Comb5="COX4",
  Comb6="COX4"
))
Combo_names_568 <- as_labeller(c(
  Comb5="gH2AX",
  Comb6="SDHB"
))
ggplot(X0Hrs_wt ,aes(x =Treatment, fill =  Treatment, y =  log10(nuclei_selected_selected_formula_nuc_cyt_488)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  facet_wrap("Combo", labeller = Combo_names_488, scales = "free_y")+ggtitle("T0 WT 488 Ratio N/C")
ggsave(here::here("Output","T0 WT 488_int.pdf"), width = 15, height = 11)
ggplot(X0Hrs_wt ,aes(x =Treatment,  fill =  Treatment,  y =  log2(nuclei_selected_selected_formula_nuc_cyt_568)))+
  
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  
  facet_wrap("Combo",labeller = Combo_names_568, scales= "free_y")+ggtitle("T0 WT 568 Ratio N/C")
ggsave(here::here("Output","T0 WT 568_int.pdf"), width = 15, height = 11)

X24Hrs_wt <- left_join(X24Hrs_wt,conditions) %>% separate(condition, sep = "_",into = c("Rep","Trep","WT","Treatment","Combo")) %>% 
  mutate(Treatment_bin = if_else(Treatment == 'DMSO',0,1)) %>% na.omit() 
COX4_wt_T0_T24_nc <- rbind(X0Hrs_wt%>% subset(Combo == "Comb5")  %>% mutate(Time = "T0hrs"),
                           X24Hrs_wt%>% subset(Combo == "Comb5")  %>% mutate(Time = "T24hrs"))
ggplot(COX4_wt_T0_T24_nc ,aes(x =Treatment,fill =  Treatment, y =  log2(nuclei_selected_selected_formula_nuc_cyt_488)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  facet_wrap("Time")+ggtitle("Cox4 T0 T24 488 Ratio")

ggsave(here::here("Output","COX4 wt T0 T24 nc_int.pdf"), width = 7, height = 9)

ggplot(COX4_wt_T0_T24_nc ,aes(x =Treatment,fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cell_alexa_488_int)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  facet_wrap("Time")+ggtitle("Cox4 T0 T24 488 cell int")

ggsave(here::here("Output","COX4 wt T0 T24 cell_int.pdf"), width = 7, height = 9)


ggplot(COX4_wt_T0_T24_nc ,aes(x =Treatment,fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_nucleus_alexa_488_int)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  facet_wrap("Time")+ggtitle("Cox4 T0 T24 488 nucl int")

ggsave(here::here("Output","COX4 wt T0 T24 nucleus_int.pdf"), width = 7, height = 9)

ggplot(X24Hrs_wt ,aes(x =Treatment,  fill =  Treatment, y =  log2(nuclei_selected_selected_formula_nuc_cyt_488)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  
  
  facet_wrap("Combo",scales = "free_y", labeller = Combo_names_488)+ggtitle("T24 WT 488 Ratio N/C")
ggsave(here::here("Output","T24 WT 488_int.pdf"), width = 15, height = 11)
ggplot(X24Hrs_wt ,aes(x =Treatment,  fill =  Treatment, y =  log2(nuclei_selected_selected_formula_nuc_cyt_568)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+

  
  facet_wrap("Combo",scales = "free_y", labeller = Combo_names_568)+ggtitle("T24 WT 568 Ratio N/C")
ggsave(here::here("Output","T24 WT 568_int.pdf"), width = 15, height = 11)

ggplot(X0Hrs_wt ,aes(x =Treatment,  fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cell_alexa_488_mean)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  
  facet_wrap("Combo", labeller = Combo_names_488, scales = "free_y")+ggtitle("T0 WT 488 Cell")
ggsave(here::here("Output","T0 WT Cell 488_int.pdf"), width = 15, height = 11)
ggplot(X0Hrs_wt ,aes(x =Treatment,  fill =  Treatment,  y =  log2(nuclei_selected_selected_intensity_cell_alexa_568_mean)))+
  
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  
  facet_wrap("Combo",labeller = Combo_names_568, scales= "free_y")+ggtitle("T0 WT 568 Cell")
ggsave(here::here("Output","T0 WT Cell 568_int.pdf"), width = 15, height = 11)

ggplot(X24Hrs_wt ,aes(x =Treatment,  fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cell_alexa_488_mean)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+
  
  
  facet_wrap("Combo",scales = "free_y", labeller = Combo_names_488)+ggtitle("T24 WT 488 Cell")
ggsave(here::here("Output","T24 WT cell 488_int.pdf"), width = 15, height = 11)
ggplot(X24Hrs_wt ,aes(x =Treatment,  fill =  Treatment, y =  log2(nuclei_selected_selected_intensity_cell_alexa_568_mean)))+
  geom_boxplot(alpha = 0.7)+ theme_bw()+theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("DMSO" =scales::muted("blue") ,
                               "Etoposide" = scales::muted("red")))+

  facet_wrap("Combo",scales = "free_y", labeller = Combo_names_568)+ggtitle("T24 WT 568 Cell")
ggsave(here::here("Output","T24 WT Cell 568_int.pdf"), width = 15, height = 11)

for(i in unique(X0Hrs_wt$Combo)){
  Combination  = Combns_wt[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
  model_0hrs_488 <- lm(nuclei_selected_selected_formula_nuc_cyt_488~Trep + Rep+ Treatment_bin , data =X0Hrs_wt %>% subset(
    Combo == i
  ) )
  # 
  summary_tbl <- summary(model_0hrs_488)
  summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
  list_of_significance[[glue::glue("X0Hrs_wt {i} {Combination[1]}")]] <- summary_tbl
  model_0hrs_568 <- lm( nuclei_selected_selected_formula_nuc_cyt_568 ~Trep +Rep+ Treatment_bin , data =X0Hrs_wt %>% subset(
    Combo == i ))
  
  summary_tbl <- summary(model_0hrs_568)
  summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
  list_of_significance[[glue::glue("X0Hrs_wt {i} {Combination[2]}")]] <- summary_tbl
  
}

for(i in na.omit(unique(X24Hrs_wt$Combo))){
  Combination  = Combns_wt[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
  model_0hrs_488 <- lm(nuclei_selected_selected_formula_nuc_cyt_488~Trep + Rep+ Treatment_bin, data =X24Hrs_wt %>% subset(
    Combo == i
  ) )
  # 
  summary_tbl <- summary(model_0hrs_488)
  summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
  list_of_significance[[glue::glue("X24Hrs_wt {i} {Combination[1]}")]] <- summary_tbl
  
  model_0hrs_568 <- lm( nuclei_selected_selected_formula_nuc_cyt_568 ~Trep +Rep+ Treatment_bin, data =X24Hrs_wt %>% subset(
    Combo == i ))
   summary_tbl <- summary(model_0hrs_568)
  summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
  list_of_significance[[glue::glue("X24Hrs_wt {i} {Combination[2]}")]] <- summary_tbl
  
}

for(i in unique(X0Hrs_wt$Combo)){
  Combination  = Combns_wt[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
  model_0hrs_488 <- lm(log2(nuclei_selected_selected_intensity_cell_alexa_488_mean)~Trep + Rep+ Treatment_bin,
                       data =X0Hrs_wt %>% subset(
                         Combo == i
                       ) )
  # 
  summary_tbl <- summary(model_0hrs_488)
  summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
  list_of_significance[[glue::glue("X0Hrs_wt Cell {i} {Combination[1]}")]] <- summary_tbl
 
  model_0hrs_568 <- lm(log2(nuclei_selected_selected_intensity_cell_alexa_568_mean) ~Trep +Rep+ Treatment_bin, data =X0Hrs_wt %>% subset(
    Combo == i ))
  
  summary_tbl <- summary(model_0hrs_568)
  summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
  list_of_significance[[glue::glue("X0Hrs_wt Cell {i} {Combination[2]}")]] <- summary_tbl
  
}

for(i in na.omit(unique(X24Hrs_wt$Combo))){
  Combination  = Combns_wt[X1 == i,2:3] %>% unlist()# %>% paste0(collapse = "_")
  model_0hrs_488 <- lm(log2(nuclei_selected_selected_intensity_cell_alexa_488_mean)~Trep + Rep+ Treatment_bin, data =X24Hrs_wt %>% subset(
    Combo == i
  ) )
  # 
  summary_tbl <- summary(model_0hrs_488)
  summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
  list_of_significance[[glue::glue("X24Hrs_wt Cell {i} {Combination[1]}")]] <- summary_tbl
  
  
  model_0hrs_568 <- lm( log2(nuclei_selected_selected_intensity_cell_alexa_568_mean) ~Trep +Rep+ Treatment_bin, data =X24Hrs_wt %>% subset(
    Combo == i ))
  summary_tbl <- summary(model_0hrs_568)
  summary_tbl <- summary_tbl$coefficients %>% as.data.frame() %>% rownames_to_column("Coef")
  list_of_significance[[glue::glue("X24Hrs_wt Cell {i} {Combination[2]}")]] <- summary_tbl
  
}
openxlsx::write.xlsx(list_of_significance,here::here("Datasets","Processed","Linear_Image Quant Models summary.xlsx"), overwrite = T)
