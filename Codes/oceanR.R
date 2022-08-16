##### oceanR
##################################### OCEAN ########################################

#'R package for metabolic enzyme enrichment analysis
#'
#'OCEAN: a method that defines metabolic enzyme footprint from a curated reduced
#'       version of the recon2 reaction network. The metabolic enzyme footprints
#'       are used to explore coordinated deregulations of metabolite abundances
#'       with respect to their position relative to metabolic enzymes.
#'       This is similar to Kinase-substrate and TF-targets enrichment analyses.

## install ocEAn
# install_github("saezlab/ocean")
pacman::p_load(devtools,
               ocean,tidyverse)



##Differential analysis
#comparisons <- list('tumorVsHealthy' = c(1,-2))
list_of_comparisons <- list(
    list('EtopvsDMSO' = c(1,-2)),
    list('T8vsT0' = c(1,-2)),
    list('T24vsT8' = c(1,-2)),
    list('T24vsT0' = c(1,-2)),
    list('T24vsDMSO' = c(1,-2)),
    list("PRDX1vsWT"= c(1,-2)),
    list("PRDX1_etopvsPRDX1_un"= c(1,-2)),
    list("PRDX1_8vsPRDX1_0"= c(1,-2)),
    list("PRDX1_24vsPRDX1_8"= c(1,-2)),
    list("PRDX1_24vsPRDX1_un" = c(1,-2))
)
list_of_comparisons <- list_of_comparisons %>% set_names(.,list_of_comparisons %>% flatten() %>% names())
Metabo_Data <-  read_csv(here::here("Datasets","Processed","Metabolomics_data_norm.csv"))%>% as.data.frame()%>% 
    column_to_rownames("...1")
list_of_Data <- list(EtopvsDMSO = Metabo_Data[,c("wt_etop_0_1","wt_etop_0_2","wt_etop_0_3","wt_ut_1","wt_ut_2","wt_ut_3")],
                     T8vsT0 = Metabo_Data[,c("wt_etop_8h_1","wt_etop_8h_2","wt_etop_8h_3","wt_etop_0_1","wt_etop_0_2","wt_etop_0_3")],
                     T24vsT8= Metabo_Data[,c("wt_etop_24h_1","wt_etop_24h_2","wt_etop_24h_3", "wt_etop_8h_1","wt_etop_8h_2","wt_etop_8h_3")],
                     T24vsT0 = Metabo_Data[,c("wt_etop_24h_1","wt_etop_24h_2","wt_etop_24h_3", "wt_etop_0_1","wt_etop_0_2","wt_etop_0_3")],
                     T24vsDMSO = Metabo_Data[,c("wt_etop_24h_1","wt_etop_24h_2","wt_etop_24h_3", "wt_ut_1","wt_ut_2","wt_ut_3")],
                     PRDX1vsWT = Metabo_Data[,c("prdx1_ko_ut_1","prdx1_ko_ut_2","prdx1_ko_ut_3","wt_ut_1","wt_ut_2","wt_ut_3")],
                     PRDX1_etopvsPRDX1_un = Metabo_Data[,c("prdx1_ko_etop_0_1","prdx1_ko_etop_0_2","prdx1_ko_etop_0_3","prdx1_ko_ut_1","prdx1_ko_ut_2","prdx1_ko_ut_3")],
                     PRDX1_8vsPRDX1_0 = Metabo_Data[,c("prdx1_ko_etop_8h_1","prdx1_ko_etop_8h_2","prdx1_ko_etop_8h_3","prdx1_ko_etop_0_1","prdx1_ko_etop_0_2","prdx1_ko_etop_0_3")],
                     PRDX1_24vsPRDX1_8 = Metabo_Data[,c( "prdx1_ko_etop_24h_1","prdx1_ko_etop_24h_2","prdx1_ko_etop_24h_3","prdx1_ko_etop_8h_1","prdx1_ko_etop_8h_2","prdx1_ko_etop_8h_3")],
                     PRDX1_24vsPRDX1_un = Metabo_Data[,c( "prdx1_ko_etop_24h_1","prdx1_ko_etop_24h_2","prdx1_ko_etop_24h_3","prdx1_ko_ut_1","prdx1_ko_ut_2","prdx1_ko_ut_3")]
                     
                     )

Ocean_Decoupler <- function(comparison, metabolite_data){
   # comparison <- list_of_comparisons$EtopvsDMSO
   # metabolite_data <- list_of_Data$EtopvsDMSO
    limmaRes <- runLimma(measurements = metabolite_data,#toy_metabolomic_data
                         targets = data.frame(sample = colnames(metabolite_data),
                                              condition =colnames(metabolite_data) %>% str_remove_all("_.$") ),#  toy_targets,
                         comparisons = comparison)
    
    ##Format differential analysis result
    t_table <- ttop_list_to_t_table(
        limma_res_to_ttop_list(limma_res = limmaRes,
                               comp_names = names(comparison),
                               number = length(metabolite_data[,1]),
                               adjust.method = "fdr"))
    
    ##This step is particularly important because this is where the users metabolic identifiers are mapped to the kegg ids used by the method.
    #Thus, the user should provide a mapping table in the same format as the mapping_table presented here (you can look at it to inspire yourself from it)
    #The mapping table should map the users own metabolic identifiers to kegg compound IDs
    t_table <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                                  mapping_table = Mapping %>% dplyr::select(-Query),
                                                  affixes = c("c","l","x","m","e","n","r"))
    
    ##Prepare the metabolic enzyme sets
    penalty_min <- 8 #minimum 1 and integer
    penalty_max <- 8 #maximum 9 and integer
    
    ##These are the available pathways to choose from
    #View(unique(recon2_redhuman$pathway))
    
    #Select pathways relevant to include in network
    TCA_network <- model_to_pathway_sif(pathway_to_keep =recon2_redhuman$pathway %>% na.omit())#  c("Citric acid cycle",
    # "Glycolysis/gluconeogenesis",
    # "Pyruvate metabolism",
    # "Transport, mitochondrial"))

    ##Translate enzyme complexes by mapping identifiers to names
    TCA_network_trans <- translate_complexes(TCA_network)

    ##This is to simplify network structure by removing cofactors
    TCA_network_nocofact <- remove_cofactors(TCA_network_trans)

    ##This is to simplify the network structure by compressing redundant transporters
    TCA_network_nocofact <- compress_transporters(sub_network_nocofact = TCA_network_nocofact)

    ##This is to avoid cross ping pong between reactants and products of reversible transaminases
    TCA_network_nocofact <- split_transaminases(sub_network_nocofact = TCA_network_nocofact)

    ##This is to filter all unique enzymes
    enzymes <- unique(TCA_network_nocofact$attributes$V1)
    enzymes <- enzymes[!grepl("_[clxmenr]$",enzymes)]
    
    #branch_length applies a cutoff on the minimum length of the reaction network
    #upstream and down stream of a given enzyme. Here we use a minimum length of 3
    #for both directions.
    # TCA_forest <- forestMaker(enzymes, TCA_network_nocofact$reaction_network,
    #                           branch_length = c(3,3))
    # saveRDS(TCA_forest, here::here("Datasets","Processed","TCA_forest.rds"))
    TCA_forest <- readRDS(here::here("Datasets","Processed","TCA_forest.rds"))
    reaction_set_list <- prepare_metabolite_set(penalty_range = penalty_min:penalty_max,  
                                                forest = TCA_forest,
                                                measured_metabolites = t_table$KEGG)
    
    reaction_set_list_merged <- condense_metabolite_set(reaction_set_list = reaction_set_list)
    
    penalty <- 8 #has to be between penalty_min and penalty_max and integer
    
    regulons_df <- prepare_regulon_df(reaction_set_list_merged, penalty, c(0,1))
    
    ##Compute metabolic enzyme enrichment score
    metactivity_res <- metactivity(metabolomic_t_table = t_table, 
                                   regulons_df = regulons_df, 
                                   compartment_pattern = "_[a-z]$", 
                                   k = 1000)
    
    mean_ES_df <- metactivity_res$ES
    mean_NES_df <- metactivity_res$NES
    
    ##translate the metabolic ids back to names
    translated_results <- translate_results(regulons_df = regulons_df,
                                            t_table = t_table,
                                            mapping_table = mapping_table)
    
    ##Visualise results for single enzymes
    # plots <- plotMetaboliteContribution(enzyme = 'PRDX1',
    #                                     stat_df = translated_results$t_table,
    #                                     metabolite_sets = translated_results$regulons_df,
    #                                     contrast_index = 1,
    #                                     stat_name = 't',
    #                                     scaling_factor = 5,
    #                                     nLabels = 15)

     # plot(plots$scatter)
     # plot(plots$cumsumPlot)
    #screen+Chromatin with metabolomics
    # Overlapping_terms= kegg_pathways %>% subset(str_detect(gene,paste(c(All_interesting_genes,Significant_genes),collapse = "|"))) %>% add_count(term)
    # write.xlsx(Overlapping_terms,here::here("Datasets","Processed","Chrom_Screen_hits_of_KEGG.xlsx"))
    #Visualise results at KEGG pathway level
    # pathway <- "KEGG_PYRIMIDINE_METABOLISM"
    # hm <- pathway_HM(mean_NES_df = mean_NES_df %>% subset(!duplicated(EtopvsDMSO)), pathway_name = pathway, pathways = kegg_pathways)
    # plot(hm)+ggtitle(glue::glue("{pathway} Etop DMSO WT"))
    # ggsave(here::here("Output","Figures",glue::glue("{pathway} Etop DMSO WT .pdf")), width = 20, height = 20)
    #Visualise the network
    # network_subset <- 
    #   mean_NES_df
    # pathway_name = pathway
    # pathways =kegg_pathways
    # sort_by = 1
    # current_pathway <- pathways[pathways$term == pathway_name, 
    #                               "gene"]
    # mean_NES_df_to_subset <- mean_NES_df
    # 
    # mean_NES_df_to_subset$ID_short <- gsub(">.*", "", mean_NES_df$KEGG)
    # mean_NES_df_to_subset$ID_short <- gsub("_.*", "", mean_NES_df_to_subset$ID_short)
    # 
    # current_pathway_activities <- mean_NES_df[mean_NES_df_to_subset$ID_short %in% 
    #                                               current_pathway, ]
    #   row.names(current_pathway_activities) <- current_pathway_activities$KEGG
    #   current_pathway_activities <- current_pathway_activities[order(current_pathway_activities[, 
    #                                                                                             sort_by + 1], decreasing = F), ]
    #   current_pathway_act_melt <- melt(current_pathway_activities)
    #   current_pathway_act_melt$KEGG <- factor(current_pathway_act_melt$KEGG, 
    #                                           levels = unique(current_pathway_act_melt$KEGG))
    #   names(current_pathway_act_melt) <- c("Enzyme", "Contrast", 
    #                                        "NES")
    # TCA_network_nocofact_pathways <- TCA_network_nocofact
    # Connected <-  TCA_network_nocofact_pathways$reaction_network  %>% subset(target %in% current_pathway_act_melt$Enzyme |source %in% current_pathway_act_melt$Enzyme) %>% unlist() %>% unique
    # TCA_network_nocofact_pathways$attributes <-TCA_network_nocofact_pathways$attributes  %>% subset(V1 %in% Connected |V2 %in% Connected)
    # TCA_network_nocofact_pathways$reaction_network <- TCA_network_nocofact_pathways$reaction_network  %>% subset(target %in% Connected |source %in% Connected)
    # 
    # plot_reaction_network(TCA_network_nocofact_pathways, t_table, mean_NES_df %>% subset(!duplicated(EtopvsDMSO)),
    #                       column_index = 1, vis.height = 2000)

    
    #### enrichment analysis ##### 
    ##Select all pathways from recon2_redhuman$pathway to be included in network
    # all_pathways <- as.vector(unique(recon2_redhuman$pathway)) %>% na.omit()
    
    ##Create a reaction network from the recon_redhuman model (metabolites & enzymes)
    # reaction_network <- model_to_pathway_sif(pathway_to_keep = all_pathways)
    
    ##Translate enzyme complexes by mapping identifiers to names
    # reaction_network <- translate_complexes(reaction_network)
    
    ##Create data frame metabolite - affiliated enzyme instead of source-target
    metabolites <- rearrange_dataframe(TCA_network)
    
    ##Remove compartment information and "cpd:" to get pure KEGG IDs
    metabolites$metabolites <- get_pure_kegg_ids(metabolites$metabolites)
    metabolites <- distinct(metabolites)  #keep only unique rows
    
    ##Map pathways to metabolites
    metabolites_pathway_df <- map_pathways_to_metabolites(metabolites)
    
    #Save metabolites-pathways data frame
    #write.csv(metabolites_pathway_df,"./results/metabolites_pathway_df.csv")
    
    
    ##Metabolite enrichment analysis with package decoupleR
    
    #'Calculate the activity enrichment score and p-value of all pathways by using
    #'the conditions in the matrix (patient samples vs metabolite expression) by
    #'calculating the mean over the expression of all metabolites.
    
    ##Prepare input data
    network <- metabolites_pathway_df    #input 1: pathways and their associated metabolites
    network$mor <- 1                     #add new column for mode of regulation
    network$likelihood <- 1              #add new column for edge likelihood
    
    t_table_kegg <- t_table
    t_table_kegg$KEGG <- get_pure_kegg_ids(t_table_kegg$KEGG) #remove compartment info
    t_table_kegg <- unique(t_table_kegg)
    row.names(t_table_kegg) <- t_table_kegg$KEGG  #convert kegg column to row names
    t_table_kegg$KEGG <- NULL            #delete kegg ids column
    t_table_kegg <- unique(t_table_kegg) #ensure only 1 kegg is mapped to 1 metabolite
    mat <- as.matrix(t_table_kegg)       #input 2: expression matrix
    
    ##Perform enrichment analysis
    enrichment_dec <- decoupleR::run_mean(mat, network,
                                          .source = .data$pathway,
                                          .target = .data$metabolites,
                                          .mor = .data$mor, 
                                          .likelihood = .data$likelihood,  
                                          times = 10000,            #number of permutations
                                          seed = 42,                #a single integer
                                          sparse = TRUE,
                                          randomize_type = "rows")  #randomize matrix
    enrichment_dec$condition <- NULL
    enrichment_dec <- as.data.frame(enrichment_dec)
    colnames(enrichment_dec) <- c("statistic", "pathway", "score", "p-value")
    metabolites_oceanR <- network[,c("pathway","metabolites")] %>% inner_join(Mapping, by = c("metabolites" = "KEGG")) %>% 
      group_by(pathway) %>% 
      add_count()
    pathw_more_2_metabo <- network[,c("pathway","metabolites")] %>% inner_join(Mapping, by = c("metabolites" = "KEGG")) %>% 
        group_by(pathway) %>% 
        add_count() %>% dplyr::select(pathway,n) %>% distinct() %>% subset(n>2) %>% pull(pathway)
    ##Keep only rows with statistic "normalized_mean"
    write.csv(metabolites_oceanR,here::here("Datasets","Processed","decoupleR_metabolites_pathways.csv"))
    enrichment_norm <- enrichment_dec[enrichment_dec$statistic != "mean", ]
    enrichment_norm <- enrichment_norm %>% subset(pathway %in% pathw_more_2_metabo)
    # enrichment_norm <- enrichment_norm %>% mutate(
    #     score = score*(-1)
    # )
    ##Visualize most significant pathways in a bar plot
    score <- 1.2  #define cutoff score
    barplot = plot_significant_pathways(enrichment_norm, score)
    
   
    ggsave(here::here("Output","Figures",glue::glue("enriched_decoupleR_pathways_{names(comparison)}.pdf")), 
           plot = barplot+ggtitle(glue::glue("Significant Pathways Between {names(comparison)}")), 
           scale = 1, dpi = 300, limitsize = TRUE)
    enrichment_norm
}
metabolites_oceanR <- read.csv(here::here("Datasets","Processed","decoupleR_metabolites_pathways.csv"))
scores <- purrr::map2(.x = list_of_comparisons, 
            .y = list_of_Data,
            ~Ocean_Decoupler(.x, .y))
scores[["PRDX1_24vsPRDX1_un"]] <- enrichment_norm
# openxlsx::write.xlsx(scores,here::here("Datasets","Processed","oceanR_scores.xlsx"),overwrite = T)
sheet_names <- openxlsx::getSheetNames(here::here("Datasets","Processed","oceanR_scores.xlsx"))
scores <- purrr::map(.x = sheet_names,
                     ~openxlsx::read.xlsx(here::here("Datasets","Processed","oceanR_scores.xlsx"), sheet = .x)) %>% 
  set_names(sheet_names)
rbind(scores$EtopvsDMSO %>% mutate(Comp = "EtopvsDMSO"),
      scores$PRDX1_etopvsPRDX1_un %>% mutate(Comp = "PRDX1EtopvsDMSO")) %>%
    dplyr::select(pathway,score,Comp) %>% 
    pivot_wider(names_from = Comp,values_from = score) %>% 
    mutate(Diff = PRDX1EtopvsDMSO -EtopvsDMSO ) %>% 
    ggplot(aes(x = EtopvsDMSO,y = PRDX1EtopvsDMSO,colour = Diff, label = pathway))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1)+
    scale_colour_gradient2(low = "blue",high = "red", mid = "grey80")+
    ggrepel::geom_label_repel(data = . %>% subset(abs(Diff)>2.3)) +
    theme_bw()+
    lims(x = c(-4,4),
         y = c(-4,4))+
  ggtitle("Metabolic Pathways Affected by Etop",
          subtitle = "Differences between PRDX1 and WT background")
ggsave(here::here("Output","Figures","OceanR_Etop_PRDX1_WT_vs_DMSO_pathways.pdf"))
rbind(scores$T24vsDMSO %>% mutate(Comp = "T24vsDMSO"),
      scores$EtopvsDMSO %>% mutate(Comp = "EtopvsDMSO")) %>%
  dplyr::select(pathway,score,Comp) %>% 
  pivot_wider(names_from = Comp,values_from = score) %>% 
  mutate(Diff = T24vsDMSO -EtopvsDMSO ) %>% 
  ggplot(aes(x = EtopvsDMSO,y = T24vsDMSO,colour = Diff, label = pathway))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  scale_colour_gradient2(low = "#053061",high = "#670A1F", mid = "grey80")+
  ggrepel::geom_text_repel(data = . %>% subset(abs(Diff)>1.3),size = 12) +
  theme_bw()+
  theme(text = element_text(size = 35))+
  lims(x = c(-3.2,3.8),
       y = c(-3.2,3.8))+
  ggtitle("Metabolic Pathways Affected by Recovery Etop",
          subtitle = "Differences between T24 and T0 background \n Pyrimidine Synthesis Seems to be related to the recovery rather than the damage")
ggsave(here::here("Output","Figures","OceanR_Etop_T24_vs_T0_pathways.pdf"),width = 15,height = 10)

rbind(scores$PRDX1_24vsPRDX1_un %>% mutate(Comp = "PRDX1_24vsPRDX1_un"),
      scores$PRDX1_etopvsPRDX1_un %>% mutate(Comp = "PRDX1_etopvsPRDX1_un")) %>%
  dplyr::select(pathway,score,Comp) %>% 
  pivot_wider(names_from = Comp,values_from = score) %>% 
  mutate(Diff = PRDX1_24vsPRDX1_un -PRDX1_etopvsPRDX1_un ) %>% 
  ggplot(aes(x = PRDX1_etopvsPRDX1_un,y = PRDX1_24vsPRDX1_un,colour = Diff, label = pathway))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  scale_colour_gradient2(low = "#053061",high = "#670A1F", mid = "grey80")+
  ggrepel::geom_text_repel(data = . %>% subset(abs(Diff)>1.3),size = 12) +
  theme_bw()+
  theme(text = element_text(size = 35))+
  
  lims(x = c(-3.2,3.8),
       y = c(-3.2,3.8))+
  ggtitle("Metabolic Pathways Affected by Recovery Etop in PRDX1 KO")
ggsave(here::here("Output","Figures","OceanR_Etop_T24_vs_T0_PRDX1_pathways.pdf"),width = 15,height = 10)

Enzyme_Importance <- mean_NES_df %>% 
    mutate(KEGG = KEGG %>% str_remove_all("_[:graph:]*$") %>% 
               str_remove_all(">[:graph:]*$")  ) %>%
    group_by(KEGG) %>% 
    summarise(Max_effect = max(abs(Etop0VsWT), na.rm = T)) %>% 
    rename(ID = KEGG) %>% 
    inner_join(Human_hsa)
KEGG_list <-Enzyme_Importance %>% 
    mutate(Max_effect = scale(Max_effect)[,1]) %>%
    arrange(-Max_effect) %>% 
    pull(Max_effect,KEGG )
kk2 <- gseKEGG(geneList     = KEGG_list,
               organism     = 'hsa',
               minGSSize    = 10,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
if((kk2@result %>% nrow)>0){
    
    enrichplot::ridgeplot(kk2, showCategory = 100)+
        ggtitle(glue::glue("OcearR enzyme important"," Gene  KEGG "))
    ggsave(here::here(output_folder, glue::glue("OcearR enzyme important"," KEGG.png")), height = 20, width  = 15)}


mkk2 <- gseMKEGG(geneList = KEGG_list,
                 organism = 'hsa',
                 minGSSize = 10,
                 pvalueCutoff = 0.05)
if((mkk2@result %>% nrow)>0){
    
    enrichplot::ridgeplot(mkk2, showCategory = 68)+
        ggtitle(glue::glue("OcearR enzyme important"," Gene MKEGG "))
    ggsave(here::here(output_folder, glue::glue("OcearR enzyme important"," MKEGG.png")), height = 20, width  = 15)}

###investigating why "Nucleotide interconversion" is positive in decoupleR Etop vs DMSO , but negative in PRDX1_EtopvsUntreated
Nucleotide_interconversion <- network %>% subset(pathway == "Nucleotide interconversion") %>% pull(metabolites)  
Nucleotide_interconversion <-     subset(Mapping, KEGG %in% Nucleotide_interconversion) %>% pull(Metabolite)
EtopvsDMSO <- list_of_Data$EtopvsDMSO[Nucleotide_interconversion,] %>% 
    rownames_to_column("Metabolites") %>% 
    pivot_longer(-Metabolites,names_to = "Condition",values_to = "Abundance") %>% 
    mutate(Condition = str_match(Condition,"^(wt_[:graph:]*?)_")[,2]) %>% 
    group_by(Condition ,Metabolites) %>% 
        summarise(Mean_ab = mean(Abundance,na.rm = T)) %>% 
    summarise(Mean_ab = mean(Mean_ab,na.rm = T)) 
PRDX_etopvsDMSO <- list_of_Data$PRDX1_etopvsPRDX1_un[Nucleotide_interconversion,] %>% 
    rownames_to_column("Metabolites") %>% 
    pivot_longer(-Metabolites,names_to = "Condition",values_to = "Abundance") %>% 
    mutate(Condition = str_match(Condition,"^(prdx1_ko_[:graph:]*?)_")[,2]) %>% 
    group_by(Condition ,Metabolites) %>% 
    summarise(Mean_ab = mean(Abundance,na.rm = T)) %>% 
    summarise(Mean_ab = mean(Mean_ab,na.rm = T)) 
rbind(EtopvsDMSO,PRDX_etopvsDMSO)
#yes! the graphs are correct. Although in wt, metabolites related to Nucleotide metabolism increase upon etop,
#in PRDX1_KO background they decrease. Which is very curious