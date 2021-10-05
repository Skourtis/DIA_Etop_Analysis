#### Etop DIA ####
pacman::p_load(piggyback, renv, here, tidyverse, targets, DEP,pheatmap,diann,PeCorA,sva,imp4p,proDA,
               org.Hs.eg.db,clusterProfiler,ggridges,SubCellBarCode,
               visNetwork,matrixStats,magick,testthat, openxlsx, janitor,seqinr)

#### Loading####
source(here::here("Codes","functions.R"))
set.seed(1234)
Proteomic_Ruler <- here::here("Datasets","Processed","CCLE_prot_Ruler.txt") %>% read.delim() %>% .[-1,] %>% 
    dplyr::select(matches("Copy|Uniprot_Acc|accuracy"))%>% 
    remove_rownames() %>% 
    column_to_rownames("Uniprot_Acc") %>% 
    #mutate(across(where(is.numeric), as.numeric)) %>% 
    set_names(.,str_remove_all(names(.),"Copy\\.number\\.")) %>% 
    mutate(across(contains("_"),~log2(as.numeric(.x))),
           across(where(is.numeric), ~if_else(is.infinite(.x), NaN,.x))) %>% 
    subset(.,rowSums(is.na(.))<(ncol(.)/3)) %>% 
    subset(Absolute.quantification.accuracy != "low") %>% 
    dplyr::select(-Absolute.quantification.accuracy) %>% 
    janitor::clean_names()

Uniprot_length_Mass <- here::here("Datasets","Raw", "Uniprot_Molecular_sizes.tab") %>% 
    read_tsv() %>% 
    dplyr::select(-`Entry name`) %>% 
    set_names(c("Uniprot","Length","Mass")) %>% 
    mutate(Mass = Mass/1000)
HUMAN_9606 <- read_tsv(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"),
                       col_names = FALSE) %>% set_names(c("Uniprot","Type","ID"))
HUMAN_9606_count <- HUMAN_9606 %>% 
    add_count(Uniprot) %>% dplyr::select(-c(Type,ID)) %>% 
    distinct()

Sabatini <- openxlsx::read.xlsx(here::here("Datasets","Raw","metabolism_gene_list_Sabatini.xlsx")) %>% 
    pull(Gene.Symbol)
Sabatini_Uniprot <- HUMAN_9606 %>% subset((ID %in% Sabatini) & Type == "Gene_Name") %>%
    subset(!duplicated(Uniprot)) %>% pull(Uniprot, ID)#%>% unique()
output_folder = here::here("Output","Etop_DIA_EC_two_methods")
Hits_Joanna <- openxlsx::read.xlsx(here::here("Datasets","Raw",
                                              "Candidates_proteins.xlsx")) %>% pull(`71.genes`) %>% c("FDX2")
Candidates_Joanna<-openxlsx::read.xlsx(here::here("Datasets","Raw",
                                                  "Candidates_proteins.xlsx"), sheet =2 , colNames = F)[,1] %>% unlist() %>%c("FDX2")

early_fast <- c("Ku80","PRKDC","XRCC4","PARP1","Ku70","LIG4",
                "PCNA","RAD50","ATM")
early_slow <- c("MDC1",'53BP1',"MRE11","BRCA1","BRCA2","RPA","BARD1","PAXIP1")
late_slow <- c("RAD54",'RAD52',"RAD51")
Proteins_known_behaviour<- data.frame(ID = c(early_fast,early_slow,late_slow),
                                      Behaviour = rep(c("early_fast","early_slow","late_slow"), c(length(early_fast),
                                                                                                  length(early_slow),
                                                                                                  length(late_slow))))
Uniprot_known_behaviour <- left_join(Proteins_known_behaviour, HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type)) %>% 
    group_split(Behaviour) %>% 
    set_names(.,map_chr(.x = .,~.x %>% pull(Behaviour) %>% unique())) %>% 
    map(~pull(.x,Uniprot) %>% unique() %>% na.omit())
KORE<- openxlsx::read.xlsx(here::here("Datasets","Raw","KORE_DDR_genes_uniprot_IDS.xlsx"))[,2]

DDR_GO<- subset(HUMAN_9606, ID %in% (read_csv(here::here("Datasets","Raw","QuickGO-annotations-DSB repair_without duplicates.csv")) %>%
                                         pull(SYMBOL))) %>% pull(Uniprot) %>% unique()
lists_of_interest = list("Early_Fast" = Uniprot_known_behaviour$early_fast,
                         "Early_Slow"=  Uniprot_known_behaviour$early_slow,
                         "Late_slow" =Uniprot_known_behaviour$late_slow,
                         "KORE" = KORE,
                         "DDR_GO" = DDR_GO,
                         "Hits_Joanna" = subset(HUMAN_9606, ID %in% (Hits_Joanna)) %>% pull(Uniprot) %>% unique(),
                         "Candidates_Joanna" = subset(HUMAN_9606, ID %in% (Candidates_Joanna)) %>% pull(Uniprot) %>% unique() )
Interesting_proteins <- data.frame(Uniprot= c(lists_of_interest %>% unlist(use.names = F))) %>% 
    mutate(Behaviour = " ")
add_multiple_behaviour <- function(list_Uniprot,behaviour){
    #list_Uniprot <- lists_of_interest[[3]]
    #behaviour <- names(lists_of_interest  )[3]
    Interesting_proteins <- Interesting_proteins %>%
        mutate(Behaviour = if_else(Uniprot %in% list_Uniprot, paste(Behaviour, behaviour, sep = ";"),Behaviour))
    Interesting_proteins
}

Interesting_proteins <- purrr::imap_dfr(lists_of_interest,add_multiple_behaviour) %>%
    distinct() %>% 
    group_by(Uniprot) %>% 
    mutate(Behaviour = paste0(Behaviour, collapse = ""),
           Behaviour = str_remove_all(Behaviour,"^ \\;")) %>% 
    distinct()

colours_compartmens <- c("#ADE5F0","#2FFEE0","#29FEFD","#34D5F7","#00D8DF","#C9B3A4","#C58F6C",
                         "#ECECEC","#D1D1D1","#A8A8A8","#6B6B6B",
                         "#FFDD17","#FFBF3D","#FFA095","#FF846D") %>% set_names(markerProteins$Compartments %>% levels()) 
contaminant_uniprots <- seqinr::read.fasta(here::here(
    "Datasets","Raw","contaminants.fasta"),seqtype ="AA") %>% names() 
Complexes <- read.delim(here::here("Datasets","Raw","coreComplexes.txt")) %>% subset(Organism == "Human") %>% 
    janitor::clean_names() %>% 
    pull(subunits_uni_prot_i_ds,complex_name)
Complexes <- Complexes %>% str_split(";") %>% set_names(names(Complexes))
#####
# DIA_report_file <- "P11833_P11841_Method_1_report.tsv"
# 
# Samples <- data.frame(run=  c("P11833","P11834","P11835","P11836","P11837" ,"P11838","P11839","P11840" ,"P11841"),
#                      
#                       
#                       Condition = c("DMSO_1","DMSO_2","DMSO_3","0H_1", "0H_2","0H_3","24H_1", "24H_2" ,"24H_3"))
# 
# 
# Method1_P11833_41 <- Load_DIA_NN_Data(DIA_report_file,Samples%>% 
#                                            mutate(across(everything(),janitor::make_clean_names)))
# Method1_P11833_41_DEP <- DEP_DIA(Method1_P11833_41,"Method1_P11833_41")
# Method1_P11833_41_DEP <- DEP_DIA(Method1_P11833_41,"Method1_P11833_41")

list_files_to_analyse <- list(DIA_report_file = "P11833_P11841_Method_1_report.tsv",
                              DIA_report_file_method2 = "P11833_P11841_Method_2_report.tsv",
                              DIA_report_file_method2_w_old = "Method2_original_report.tsv",
                              DIA_report_file_P11685 = "P11685_Discovery_report.tsv")

Samples <- data.frame(run=  c("P11833","P11834","P11835","P11836","P11837" ,"P11838","P11839","P11840" ,"P11841",
                              "p11591","p11592","p11593","p11685","p11688" ,"p11691","p11686","p11689" ,"p11692"),
                      Condition = c("DMSO_1","DMSO_2","DMSO_3","0H_1", "0H_2","0H_3","24H_1", "24H_2" ,"24H_3",
                                    "DMSO_4","DMSO_5","DMSO_6","0H_4", "0H_5","0H_6","24H_4", "24H_5" ,"24H_6"))%>% 
    mutate(across(everything(),janitor::make_clean_names))

DIA_temp <- map(.x = list_files_to_analyse[[3]], 
                ~Load_DIA_NN_Data(.x,Samples))
#It's interesting to see that a lot of swissprot uniprots not captured by CCLE are nuclear proteins, can we show statistical significance?
DIA_temp_temp <- DEP_DIA(DIA_temp,"Temporary")
Temp_proteomic_ruler_sample <- input_matrix %>%
    as.data.frame() %>% 
    rownames_to_column("PG") %>%
    mutate(Uniprot = str_remove_all(PG, ";[:graph:]*$"))
Temp_proteomic_ruler_sample_2 <- Temp_proteomic_ruler_sample %>% left_join(HUMAN_9606 %>% subset(Type == "Gene_Name"), by =  "Uniprot")%>%
    left_join(markerProteins[,1:2], by = c("ID" = "Proteins")) %>% 
    pivot_longer(cols = contains("_"), names_to = "Condition", values_to = "Abundance")
    
Temp_proteomic_ruler <-     left_join(Temp_proteomic_ruler_sample,Proteomic_Ruler %>% rownames_to_column("Uniprot") %>% 
    dplyr::select(Uniprot,u2os_bone)) %>% 
    mutate(across(contains("_"), ~.x-u2os_bone)) 
Missing_U2OS_test <- list(Missing_U2OS = Temp_proteomic_ruler %>% subset(is.na(u2os_bone)) %>% pull(Uniprot) %>% unique(),
                     Universe = Temp_proteomic_ruler %>% pull(Uniprot) %>% unique())

ego <- enrichGO(gene          = Missing_U2OS_test$Missing_U2OS,
                universe      = Missing_U2OS_test$Universe,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                keyType = "UNIPROT",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)
#enrichplot::upsetplot(ego %>% simplify(), n  = 20)
dotplot(ego, showCategory=30) + ggtitle("")

                     
Temp_proteomic_ruler <- na.omit(Temp_proteomic_ruler) %>% 
    left_join(HUMAN_9606 %>% subset(Type == "Gene_Name"), by =  "Uniprot")%>%
    left_join(markerProteins[,1:2], by = c("ID" = "Proteins")) %>% 
    pivot_longer(cols = contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
    subset(Condition != "u2os_bone")
ggplot(Temp_proteomic_ruler %>% subset(Condition == "dmso_1"),
       aes(x = reorder(PG, Abundance), y= Abundance, colour = Compartments))+
    geom_col()+
    facet_wrap("Compartments")+
    ggtitle("Perc_on chrom DMSO_1")
dmso_gsea <- gseGO(geneList     = Temp_proteomic_ruler %>% subset(Condition == "dmso_1" & !is.na(Abundance) & !duplicated(Uniprot)) %>% arrange(-Abundance) %>%
                      pull(Abundance,Uniprot) %>%  scale() %>% .[,1] ,
      OrgDb        = org.Hs.eg.db,
      ont          = "CC",
      minGSSize    = 100,
      keyType = "UNIPROT",
      maxGSSize    = 500,
      pvalueCutoff = 0.05,
      verbose      = FALSE)
ridgeplot(dmso_gsea)+ggtitle("enrichment of GO_CC in dmso_1 normalised to U2OS cells ")
#enrichplot::upsetplot(ego %>% simplify(), n  = 20)
dmso_gsea_sample <- gseGO(geneList     = Temp_proteomic_ruler_sample_2 %>% subset( Condition == "dmso_1" & !is.na(Abundance) & !duplicated(Uniprot)) %>% arrange(-Abundance) %>%
                       pull(Abundance,Uniprot) %>%  scale() %>% .[,1] ,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "CC",
                   minGSSize    = 100,
                   keyType = "UNIPROT",
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
ridgeplot(dmso_gsea_sample)+ggtitle("enrichment of GO_CC in dmso_1 unnormalised to U2OS cells ")


Temp_proteomic_ruler %>%  
    group_by(Uniprot,ID,Compartments ) %>%
    dplyr::summarise(Mean_Abundance_u2os = mean(Abundance, na.rm = T),.groups = "keep") %>% 
    ungroup() %>% 
    left_join(Temp_proteomic_ruler_sample_2 %>% 
                  group_by(Uniprot,ID,Compartments ) %>%
                  dplyr::summarise(Mean_Abundance_sample = mean(Abundance, na.rm = T),.groups = "keep") %>% 
                  ungroup()) %>% 
    mutate(across(contains("Abundance"), scale)) %>% 
    # subset(is.na(Compartments)) %>% 
    ggplot(aes(x = Mean_Abundance_u2os, y = Mean_Abundance_sample, label = ID))+
    geom_point()+
    ggrepel::geom_text_repel(data = . %>% subset(str_detect(ID, "^RPL")),max.overlaps = 50 )+
    geom_abline(slope=1, intercept = 0)+
    facet_wrap("Compartments")
    

    
Output_proteomic_ruler(data_norm@assays@data@listData[[1]] %>% as.data.frame(), "testing_temp", Uniprot_length_Mass )
DIA_matrices <- map(.x = list_files_to_analyse, 
                    ~Load_DIA_NN_Data(.x,Samples%>% 
                                          mutate(across(everything(),janitor::make_clean_names))))
DIA_DEP <- map2(DIA_matrices,list_files_to_analyse,DEP_DIA)


P11685$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot", 'significant')] %>%
    dplyr::rename(P11685_Log2_FC_24_vs_0 = log2_FC,
                  significant_P11685 = significant) %>% 
    inner_join(Method2_P11833_41_DEP$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot", 'significant')] %>% 
                   dplyr::rename(Method2_Log2_FC_24_vs_0 = log2_FC,
                                 significant_method2 = significant)) %>% 
    mutate(significant = case_when(
        significant_method2 == T  & significant_P11685==T ~ "Both",
        significant_method2 == T  & significant_P11685==F ~ "method2",
        significant_method2 == F  & significant_P11685==T ~ "P11685",
        significant_method2 == F  & significant_P11685==F ~ "neither"
        
    )) %>% 
    ggplot(aes(x = P11685_Log2_FC_24_vs_0, y =  Method2_Log2_FC_24_vs_0,colour =  significant ))+
    geom_point(alpha = 0.1, data = . %>% subset(significant == "neither"))+
    geom_point(alpha = 1, data = . %>% subset(significant != "neither"))+
    theme_bw()+
    ggtitle("Correlation of Fold change between P11685 Samples and method 2 P11833 samples",
            subtitle = "24hrs vs 0 hrs chromatin")
    
Method1_P11833_41_DEP$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot")] %>%
    dplyr::rename(Method1_Log2_FC_24_vs_0 = log2_FC) %>% 
    inner_join(Method2_P11833_41_DEP$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot")] %>% 
                   dplyr::rename(Method2_Log2_FC_24_vs_0 = log2_FC)) %>% 
    ggplot(aes(x = P11685_Log2_FC_24_vs_0, y =  Method2_Log2_FC_24_vs_0))+
    geom_point()

edata = cbind(DIA_DEP$DIA_report_file$ImputtedMethod2_P11833_41_DEP$Unimputted %>% na.omit())
batch = if_else(str_match(colnames(Method2_P11833_41_DEP$Imputted),"_(.)$")[,2] %>% as.numeric() <4,2,1)


# parametric adjustment
combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

pca_res <- prcomp(combat_edata1  %>% na.omit() %>% t(), scale=TRUE) 
#plot_missval(data_filt)
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)

var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)

pca_res$x %>% 
    as.data.frame %>%
    rownames_to_column("Sample") %>% 
    mutate(Condition = str_remove(Sample,"_.$")) %>% 
    ggplot(aes(x=PC1,y=PC2, label = Sample, colour = Condition )) + geom_point(size=4) +
    ggrepel::geom_label_repel()+
    labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
         y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
    theme(legend.position="top") +
    ggtitle(dataset_name)+ 
    theme(plot.title = element_text(size = 20))

Pecora_result <- PeCoRa_function(list_files_to_analyse$DIA_report_file_method2, Samples)

PeCorA::PeCorA_plotting(Pecora_result$disagree_peptides,Pecora_result$disagree_peptides[1,],Pecora_result$scaled_peptides)
row.names(testing) <- testing$peptide
ego <- enrichGO(gene          = testing %>% subset(adj_pval<0.01) %>% pull(protein) %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-.") %>% unique(),
                universe      = testing %>% pull(protein) %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-.")%>% unique(),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                keyType = "UNIPROT",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
#enrichplot::upsetplot(ego %>% simplify(), n  = 20)
dotplot(ego%>% simplify(), showCategory=30) + ggtitle("dotplot for ORA")
hist(log2(df_input[, "Precursor.Normalised"]), 100, main = "Histogram of log2 intensities", 
     col = "steelblue", border = "steelblue", freq = FALSE)





