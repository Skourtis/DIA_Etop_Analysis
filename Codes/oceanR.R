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
               ocean)



##Differential analysis
#comparisons <- list('tumorVsHealthy' = c(1,-2))
comparisons <- list('Etop0VsWT' = c(1,-2))
metabolomics_wt_etop <- read_csv(here::here("Datasets","Processed","Metabolomics_data_norm.csv"))%>% as.data.frame()%>% 
    column_to_rownames("...1") %>% 
    dplyr::select(!matches("prdx|noco|8|24"))
limmaRes <- runLimma(measurements = metabolomics_wt_etop,#toy_metabolomic_data
                     targets = data.frame(sample = colnames(metabolomics_wt_etop),
                                          condition =colnames(metabolomics_wt_etop) %>% str_remove_all("_.$") ),  #toy_targets,
                     comparisons = comparisons)

##Format differential analysis result
t_table <- ttop_list_to_t_table(
    limma_res_to_ttop_list(limma_res = limmaRes,
                           comp_names = names(comparisons),
                           number = length(metabolomics_wt_etop[,1]),
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
View(unique(recon2_redhuman$pathway))

##Select pathways relevant to include in network
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
TCA_forest <- forestMaker(enzymes, TCA_network_nocofact$reaction_network,
                          branch_length = c(3,3))
saveRDS(TCA_forest, here::here("Datasets","Processed","TCA_forest.rds"))
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
plots <- plotMetaboliteContribution(enzyme = 'ALDH4A1', 
                                    stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_results$regulons_df, 
                                    contrast_index = 1, 
                                    stat_name = 't', 
                                    scaling_factor = 5,
                                    nLabels = 15)

plot(plots$scatter)
plot(plots$cumsumPlot)

##Visualise results at KEGG pathway level
hm <- pathway_HM(mean_NES_df = mean_NES_df, pathway_name = 'KEGG_PROPANOATE_METABOLISM', pathways = kegg_pathways)
plot(hm)

##Visualise the network
plot_reaction_network(TCA_network_nocofact, t_table, mean_NES_df,
                      column_index = 1, vis.height = 2000)


#### enrichment analysis ##### 
##Select all pathways from recon2_redhuman$pathway to be included in network
all_pathways <- as.vector(unique(recon2_redhuman$pathway)) %>% na.omit()

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
enrichment <- decoupleR::run_mean(mat, network,
                       .source = .data$pathway,
                       .target = .data$metabolites,
                       .mor = .data$mor, 
                       .likelihood = .data$likelihood,  
                       times = 10000,            #number of permutations
                       seed = 42,                #a single integer
                       sparse = TRUE,
                       randomize_type = "rows")  #randomize matrix
enrichment$condition <- NULL
enrichment <- as.data.frame(enrichment)
colnames(enrichment) <- c("statistic", "pathway", "score", "p-value")

##Keep only rows with statistic "normalized_mean"
enrichment_norm <- enrichment[enrichment$statistic != "mean", ]

##Visualize most significant pathways in a bar plot
score <- 0.8  #define cutoff score
barplot = plot_significant_pathways(enrichment_norm, score)
barplot

ggsave(here::here("Output","Metabolomics","enriched_decoupleR_pathways.png"), plot = barplot, scale = 1, dpi = 300, limitsize = TRUE)

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

list(Enzymes_oceanR = Enzyme_Importance$ID,
     Chromatin_Etop = Methods_DIA$)