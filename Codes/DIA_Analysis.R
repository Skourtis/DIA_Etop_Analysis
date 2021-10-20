#### Etop DIA ####
pacman::p_load(piggyback, renv, here, tidyverse, targets, DEP,pheatmap,diann,PeCorA,sva,imp4p,proDA,
               org.Hs.eg.db,clusterProfiler,ggridges,SubCellBarCode,eulerr,
               visNetwork,matrixStats,magick,testthat, openxlsx, janitor,seqinr)
#Change the inference to Protein level
#Double run through NN
#### Loading####
source(here::here("Codes","functions.R"))
Human_GEM_genes <- here::here("Datasets","Raw","Human_GEM_model_genes.tsv") %>% read_tsv() %>% 
  pull(geneSymbols) %>% unique()
set.seed(1234)
Proteomic_Ruler <- here::here("Datasets","Processed","CCLE_prot_Ruler.txt") %>% read.delim() %>% .[-1,] %>% 
    dplyr::select(matches("Copy|Uniprot_Acc|accuracy"))%>% 
    remove_rownames() %>% 
    column_to_rownames("Uniprot_Acc") %>% 
    #mutate(across(where(is.numeric), as.numeric)) %>% 
    set_names(.,str_remove_all(names(.),"Copy\\.number\\.")) %>% 
    mutate(across(contains("_"),~log2(as.numeric(.x))),
           across(where(is.numeric), ~if_else(is.infinite(.x), NaN,.x))) %>% 
    subset(!is.nan(U2OS_BONE)) %>% 
    # subset(.,rowSums(is.na(.))<(ncol(.)/3)) %>%
     subset(Absolute.quantification.accuracy != "low") %>%
     dplyr::select(-Absolute.quantification.accuracy) %>%
    janitor::clean_names()
# Proteomic_Ruler_acc <- here::here("Datasets","Processed","CCLE_prot_Ruler.txt") %>% read.delim() %>% .[-1,] %>% 
#     dplyr::select(matches("Copy|Uniprot_Acc|accuracy"))%>% 
#     remove_rownames() %>% 
#     column_to_rownames("Uniprot_Acc") %>% 
#     #mutate(across(where(is.numeric), as.numeric)) %>% 
#     set_names(.,str_remove_all(names(.),"Copy\\.number\\.")) %>% 
#     mutate(across(contains("_"),~log2(as.numeric(.x))),
#            across(where(is.numeric), ~if_else(is.infinite(.x), NaN,.x))) %>% 
#     subset(.,rowSums(is.na(.))<(ncol(.)/3)) %>% 
#     subset(Absolute.quantification.accuracy != "low") %>% 
#     # dplyr::select(-Absolute.quantification.accuracy) %>% 
#     janitor::clean_names()

Uniprot_length_Mass <- here::here("Datasets","Raw", "Uniprot_Molecular_sizes.tab") %>% 
    read_tsv() %>% 
    dplyr::select(-`Entry name`) %>% 
    set_names(c("Uniprot","Length","Mass")) %>% 
    mutate(Mass = Mass/1000)
HUMAN_9606 <- read_tsv(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"),
                       col_names = FALSE) %>% set_names(c("Uniprot","Type","ID"))
Human_hsa <- HUMAN_9606 %>% 
  subset(Type == "KEGG") %>% 
  mutate(ID = str_remove_all(ID,"hsa:")) %>% 
  dplyr::select(-Type) %>% 
  subset(!duplicated(Uniprot))
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
list_of_metabolism = list("Metabolic_library" = Sabatini_Uniprot,
                          "Human_GEM" = HUMAN_9606 %>% 
                            subset(ID %in% Human_GEM_genes) %>% 
                            pull(Uniprot) %>% unique())

lists_of_interest = list("Early_Fast" = Uniprot_known_behaviour$early_fast,
                         "Early_Slow"=  Uniprot_known_behaviour$early_slow,
                         "Late_slow" =Uniprot_known_behaviour$late_slow,
                         "KORE" = KORE,
                         "DDR_GO" = DDR_GO,
                         "Hits_Joanna" = subset(HUMAN_9606, ID %in% (Hits_Joanna)) %>% pull(Uniprot) %>% unique(),
                         "Candidates_Joanna" = subset(HUMAN_9606, ID %in% (Candidates_Joanna)) %>% pull(Uniprot) %>% unique() )

add_multiple_behaviour <- function(list_Uniprot,behaviour){
    # list_Uniprot <- list_of_metabolism[[1]]
    # behaviour <- names(list_of_metabolism  )[1]
    Interesting_proteins <- data.frame(Uniprot= c(list_Uniprot %>% unlist(use.names = F))) %>% 
      mutate(Behaviour = " ") %>%
        mutate(Behaviour = if_else(Uniprot %in% list_Uniprot, paste(Behaviour, behaviour, sep = ";"),Behaviour))
    Interesting_proteins
}
Metabolic_proteins <- purrr::imap_dfr(list_of_metabolism,add_multiple_behaviour) %>%
  distinct() %>% 
  group_by(Uniprot) %>% 
  mutate(Behaviour = paste0(Behaviour, collapse = ""),
         Behaviour = str_remove_all(Behaviour,"^ \\;")) %>% 
  distinct()

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
list_files_to_analyse <- list(DIA_report_file_1 = "P11833_P11841_Method_1_report.tsv",
                              DIA_report_file_method2 = "P11833_P11841_Method_2_report.tsv",
                              DIA_report_file_method2_w_old = "Method2_original_report.tsv",
                              DIA_report_file_method2_old_better = "Method2_original_report_new.tsv",
                              DIA_report_file_P11685 = "P11685_Discovery_report.tsv")

Samples <- data.frame(run=  c("P11833","P11834","P11835","P11836","P11837" ,"P11838","P11839","P11840" ,"P11841",
                              "p11591","p11592","p11593","p11685","p11688" ,"p11691","p11686","p11689" ,"p11692"),
                      Condition = c("DMSO_1","DMSO_2","DMSO_3","0H_1", "0H_2","0H_3","24H_1", "24H_2" ,"24H_3",
                                    "DMSO_4","DMSO_5","DMSO_6","0H_4", "0H_5","0H_6","24H_4", "24H_5" ,"24H_6"))%>% 
    dplyr::mutate(dplyr::across(tidyselect::everything(),janitor::make_clean_names))

# Combined_batch <- map(.x = list_files_to_analyse[[3]], 
#                  ~Load_DIA_NN_Data(.x,Samples))
Methods_DIA <- map(.x = list_files_to_analyse[4],
                      ~Load_DIA_NN_Data(.x,Samples))
# Methods_1_report <- read_tsv(here::here("Datasets", "Raw",list_files_to_analyse$DIA_report_file_1))
# Methods_1_report %>% 
#   janitor::clean_names() %>% 
#   subset(str_detect(run,"P11833")) %>% 
#   group_by(protein_names) %>% 
#   glimpse() %>% 
#   summarise(Amount_of_evidence = n()) %>% 
#   mutate(Shared = if_else(str_detect(protein_names,";"),"Shared","Unique")) %>% 
#  # View()
#   ggplot(aes(x = Shared, y = log2(Amount_of_evidence), colour = Shared, label =protein_names )) +
#   geom_boxplot()+
#   geom_point()+
#   ggrepel::geom_label_repel(data = . %>% subset(log2(Amount_of_evidence)>7.5))
# Methods_DIA$DIA_report_file_1 %>% row.names() %>% str_detect(";") %>% table()
# plot(eulerr::euler(list("Unique_groups" =Methods_DIA$DIA_report_file_1 %>% row.names() %>% str_subset(";", negate = T),
#                          "Shared_groups" = Methods_DIA$DIA_report_file_1 %>% row.names() %>% str_subset(";") %>% 
#                            str_split(";",simplify = T) %>% .[.%in% (Methods_DIA$DIA_report_file_1 %>% row.names() %>% str_subset(";", negate = T))]  %>% unique()),
#                     shape = "ellipse"), quantities = TRUE, main = "Uniprots in Unique and Shared Protein entries from DIANN")
# 

Methods_DIA_DEP <- map2(.x = Methods_DIA ,
                        .y = list_files_to_analyse[4],
                   ~DEP_DIA(.x,.y))
Methods_DIA_DEP$DIA_report_file_method2_old_better[-3] %>% write.xlsx(overwrite = T, file = here::here("Datasets","Processed","Batch_corrected_matrices.xlsx"), rowNames = T)
Methods_DIA_DEP$DIA_report_file_method2_old_better$DEPs %>% write.xlsx(overwrite = T, file = here::here("Datasets","Processed","Batch_corrected_DEPs.xlsx"),rowNames = T)

#testing <- Load_DIA_NN_Data(list_files_to_analyse$DIA_report_file_method2_w_old,Samples)
pept_condition <- map(.x = list_files_to_analyse, ~
                        read_tsv(here::here("Datasets","Processed",glue::glue("peptide_condition",.x,".txt")))) 
Proteomic_Ruler_input <- Import_proteomic_ruler("From_proteomic_ruler.txt")
Temp_proteomic_ruler_sample <-Proteomic_Ruler_input %>%
    as.data.frame() %>% 
    rownames_to_column("PG") %>%
    mutate(Uniprot = str_remove_all(PG, ";[:graph:]*$"))
Temp_proteomic_ruler_sample_2 <- Temp_proteomic_ruler_sample %>% 
    left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
                  subset(!duplicated(Uniprot)), by =  "Uniprot")%>%
    left_join(markerProteins[,1:2], by = c("ID" = "Proteins")) %>% 
    pivot_longer(cols = contains("_"), names_to = "Condition", values_to = "Abundance")
ggplot(Temp_proteomic_ruler_sample_2 %>% subset(Condition == "dmso_1"),
       aes(x = reorder(PG, Abundance), y= Abundance, colour = Compartments))+
    geom_col()+
    facet_wrap("Compartments")+
    ggtitle("Copy Numbers on chrom DMSO_1")

Temp_proteomic_ruler <-     left_join(Temp_proteomic_ruler_sample,Proteomic_Ruler %>% rownames_to_column("Uniprot") %>% 
    dplyr::select(Uniprot,u2os_bone)) %>% 
    mutate(across(where(is.numeric), ~.x-u2os_bone)) %>% 
  subset(!is.na(u2os_bone))
DMSO_ruler <- Temp_proteomic_ruler %>% 
  left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
              subset(!duplicated(Uniprot)), by =  "Uniprot") %>% 
 
  dplyr::select(matches("dmso|ID"))
DMSO_ruler$Mean_DMSO = DMSO_ruler %>% dplyr::select(-ID) %>%as.matrix() %>%   rowMeans2 (.,na.rm = T) %>% scale()
DMSO_ruler <- dplyr::select(DMSO_ruler, ID, Mean_DMSO)
# Missing_U2OS_test <- list(Missing_U2OS = Temp_proteomic_ruler %>% subset(is.na(u2os_bone)) %>% pull(Uniprot) %>% unique(),
#                           # Missing_Etop_DIA =  Temp_proteomic_ruler %>% subset(!is.na(u2os_bone)) %>% pull(Uniprot) %>% unique()
#                      Universe = Temp_proteomic_ruler %>% pull(Uniprot) %>% unique())
# 
# ego <- enrichGO(gene          = Missing_U2OS_test$Missing_U2OS,
#                 universe      = Missing_U2OS_test$Universe,
#                 OrgDb         = org.Hs.eg.db,
#                 ont           = "CC",
#                 keyType = "UNIPROT",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.05,
#                 qvalueCutoff  = 0.05,
#                 readable      = TRUE)
# #enrichplot::upsetplot(ego %>% simplify(), n  = 20)
# dotplot(ego, showCategory=30) + ggtitle("Proteins Found by Chromatin DIA but not U2OS CCLE_TMT")
Proteomic_Ruler %>% 
    dplyr::select(absolute_quantification_accuracy, u2os_bone ) %>% 
    rownames_to_column("Uniprot") %>% 
    na.omit() %>% 
    left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
                  subset(!duplicated(Uniprot)), by =  "Uniprot")%>%
    left_join(markerProteins[,1:2], by = c("ID" = "Proteins")) %>% 
    # pivot_longer(cols = where(is.numeric), names_to = "Condition", values_to = "Abundance") %>% 
    subset(absolute_quantification_accuracy != "low") %>% 
        ggplot( aes(x = reorder(Uniprot, u2os_bone), y= u2os_bone, colour = Compartments))+
        geom_col()+
        facet_wrap("Compartments")+
        ggtitle("WCE copy numbers u2os")
    
                     
Temp_proteomic_ruler <- Temp_proteomic_ruler %>% 
  subset(!is.na(u2os_bone)) %>% 
    left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% 
                  subset(!duplicated(Uniprot)), by =  "Uniprot")%>%
    left_join(markerProteins[,1:2], by = c("ID" = "Proteins")) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Condition", values_to = "Abundance") %>% 
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
    ggplot(aes(x = Mean_Abundance_u2os, y = Mean_Abundance_sample, label = ID, colour = Compartments))+
    geom_point()+
    ggrepel::geom_text_repel(data = . %>% subset(str_detect(ID, "^RPL")),max.overlaps = 50 )+
    geom_abline(slope=1, intercept = 0)+
    facet_wrap("Compartments")+
    ggtitle("% of Protein on chromatin vs Copy numbers on chromatin",
            "Using proteomic Ruler normalised to U2OS CCLE_TMT")
Temp_proteomic_ruler %>%  
    group_by(Uniprot,ID,Compartments ) %>%
    dplyr::summarise(Mean_Abundance_u2os = mean(Abundance, na.rm = T),.groups = "keep") %>% 
    ungroup() %>% 
    left_join(Temp_proteomic_ruler_sample_2 %>% 
                  group_by(Uniprot,ID,Compartments ) %>%
                  dplyr::summarise(Mean_Abundance_sample = mean(Abundance, na.rm = T),.groups = "keep") %>% 
                  ungroup()) %>% 
    mutate(across(contains("Abundance"), scale)) %>% 
    subset(Compartments == "N1") %>% 
    ggplot(aes(x = Mean_Abundance_u2os, y = Mean_Abundance_sample, label = ID))+
    geom_point()+
    ggrepel::geom_text_repel(data = . %>% subset(str_detect(ID, "^RPL")),max.overlaps = 50 )+
    geom_abline(slope=1, intercept = 0)+
    facet_wrap("Compartments")+
    ggtitle("% of Protein on chromatin vs Copy numbers on chromatin",
            "Using proteomic Ruler normalised to U2OS CCLE_TMT")

    
Output_proteomic_ruler(data_norm@assays@data@listData[[1]] %>% as.data.frame(), "testing_temp", Uniprot_length_Mass )
DIA_matrices <- map(.x = list_files_to_analyse, 
                    ~Load_DIA_NN_Data(.x,Samples%>% 
                                          mutate(across(everything(),janitor::make_clean_names))))
DIA_DEP <- map2(DIA_matrices,list_files_to_analyse,DEP_DIA)


Methods_DIA_DEP$DIA_report_file_P11685$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot", 'significant')] %>%
    dplyr::rename(P11685_Log2_FC_24_vs_0 = log2_FC,
                  significant_P11685 = significant) %>% 
    inner_join( Methods_DIA_DEP$DIA_report_file_method2$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot", 'significant')] %>% 
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
  geom_point(alpha = 1, size = 30, data = . %>% subset(significant == "both"))+
  
    theme_bw()+
    ggtitle("Correlation of Fold change between P11685 Samples and method 2 P11833 samples",
            subtitle = "24hrs vs 0 hrs chromatin")
    
Method1_P11833_41_DEP$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot")] %>%
    dplyr::rename(Method1_Log2_FC_24_vs_0 = log2_FC) %>% 
    inner_join(Method2_P11833_41_DEP$DEPs$x0h_vs_x24h_diff[,c("log2_FC","Uniprot")] %>% 
                   dplyr::rename(Method2_Log2_FC_24_vs_0 = log2_FC)) %>% 
    ggplot(aes(x = P11685_Log2_FC_24_vs_0, y =  Method2_Log2_FC_24_vs_0))+
    geom_point()


Pecora_result <- PeCoRa_function(list_files_to_analyse$DIA_report_file_method2, Samples)

ego <- enrichGO(gene          = Pecora_result$disagree_peptides %>% subset(adj_pval<0.05) %>% pull(protein) %>% 
                  str_remove_all(";[:graph:]*$"),
                universe      = Pecora_result$disagree_peptides  %>% pull(protein) %>% 
                  str_remove_all(";[:graph:]*$"),
                OrgDb         = org.Hs.eg.db,
                keyType = "UNIPROT",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = FALSE)
barplot(ego, showCategory = 100)+ggtitle("Protein Groups with Significant PeCoRA Peptides")

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



plot(eulerr::euler(list("Human_GEM" = Human_GEM_genes,
                  "Metabolic_Screen" = Sabatini),
             shape = "ellipse"), quantities = TRUE)
ego <- enrichGO(gene          = setdiff(Human_GEM_genes,Sabatini),
                universe      = unique(c(Human_GEM_genes,Sabatini)),
                OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                keyType = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = FALSE)
barplot( clusterProfiler::simplify(ego),showCategory=20) 

