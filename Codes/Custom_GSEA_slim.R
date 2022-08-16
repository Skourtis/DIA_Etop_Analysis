#slim Go terms#
#BiocManager::install("msigdb",force =T)
pacman::p_load("biomaRt","msigdb")

GO_annotations <- read.xlsx(here::here("Datasets","Raw","Slim_GO.xlsx"))
mart <- try(useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"))
curlHandle <- RCurl::getCurlHandle()
filters <- listFilters(mart)
attrib <-  listAttributes(mart)
#"TERM2GENE is a data.frame with first column of term ID and second column of corresponding mapped 
#gene and TERM2NAME is a data.frame with first column of term ID and second column of corresponding term name. TERM2NAME is optional."
#
# results_2 <-purrr::map_dfr(.x = ,

Go_slim <-                          getBM(attributes = c("uniprotswissprot","goslim_goa_accession", "goslim_goa_description" ) , filters = "goslim_goa",
                                  values = GO_annotations$goslim_goa_accession , mart = mart, curl=curlHandle) %>% mutate(uniprotswissprot = as.character(uniprotswissprot),
                                                                                         goslim_goa_accession = as.character(goslim_goa_accession),
                                                                                         goslim_goa_description = as.character(goslim_goa_description))
results <- Go_slim %>% subset(uniprotswissprot != "")
TERM2GENE_Slim <- results %>% dplyr::select(goslim_goa_accession,uniprotswissprot ) %>% distinct() %>% set_names(c("ID","Uniprot")) %>% 
    left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% subset(!duplicated(Uniprot)) %>% dplyr::select(-Type), by = c("Uniprot")) %>% 
    dplyr::select(-Uniprot) %>% distinct() %>% set_names(c("ID","Gene"))
TERM2NAME <- results %>% dplyr::select(goslim_goa_accession,goslim_goa_description ) %>% distinct() %>% set_names(c("ID","Name"))
em2 <- GSEA(Volcano_DFs$dmso_vs_x0h_diff %>% pull(log2_FC,ID) %>% sort(decreasing = T) ,
            TERM2GENE = TERM2GENE_Slim,
            TERM2NAME = TERM2NAME)
Slim_enchrichment_all <- function(df_input){
    GSEA(df_input %>% pull(log2_FC,ID) %>% sort(decreasing = T) ,
         TERM2GENE = TERM2GENE_Slim,
         TERM2NAME = TERM2NAME)
}
Volcano_DFs[["D14_BP"]] <- norm_genes %>% dplyr::rename(ID = Gene,
                                                        log2_FC =Etop_vs_DMS0 )
Volcano_DFs[["D10_high_BP"]] <- Screens_new$Etop_High_low_10_cnv_negsg.mle.gene_summary.txt %>% dplyr::rename(ID = Gene,
                                                                                                         log2_FC =`high|beta` )

enrichment <- map(Volcano_DFs,Slim_enchrichment_all)
for(i in 1:length(enrichment[-5])){
    enrichplot::dotplot(enrichment[[i]], showCategory  = 20)+ ggtitle(names(enrichment)[i])
}
