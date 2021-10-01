#### Etop DIA ####
#### Loading####
source(here::here("Codes","functions.R"))
set.seed(1234)
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
DIA_report_file <- "P11833_P11841_Method_1_report.tsv"

Samples <- data.frame(run=  c("P11833","P11834","P11835","P11836","P11837" ,"P11838","P11839","P11840" ,"P11841"),
                     
                      
                      Condition = c("DMSO_1","DMSO_2","DMSO_3","0H_1", "0H_2","0H_3","24H_1", "24H_2" ,"24H_3"))


Method1_P11833_41 <- Load_DIA_NN_Data(DIA_report_file,Samples%>% 
                
                                                                     mutate(across(everything(),janitor::make_clean_names)))

list_files_to_analyse <- list(DIA_report_file = "P11833_P11841_Method_1_report.tsv",
                              DIA_report_file_method2 = "P11833_P11841_Method_2_report.tsv",
                              DIA_report_file_method2_w_old = "Method2_and_original_Etop_report.tsv",
                              DIA_report_file_P11685 = "P11685_Discovery_report.tsv")

Samples <- data.frame(run=  c("P11833","P11834","P11835","P11836","P11837" ,"P11838","P11839","P11840" ,"P11841",
                              "p11591","p11592","p11593","p11685","p11688" ,"p11691","p11686","p11689" ,"p11692"),
                      Condition = c("DMSO_1","DMSO_2","DMSO_3","0H_1", "0H_2","0H_3","24H_1", "24H_2" ,"24H_3",
                                    "DMSO_4","DMSO_5","DMSO_6","0H_4", "0H_5","0H_6","24H_4", "24H_5" ,"24H_6"))

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

edata = Method2_P11833_41_DEP$Unimputted %>% na.omit()
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

Pecora_result <- PeCoRa_function(DIA_report_file, Samples)

PeCorA::PeCorA_plotting(Pecora_result$disagree_peptides,Pecora_result$disagree_peptides[,],Pecora_result$scaled_peptides)
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





