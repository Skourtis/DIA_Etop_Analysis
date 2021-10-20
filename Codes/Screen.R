#Screen
list_of_files <- c("14d-etop_vs_14d_lfc(median)-pvalue.txt",
                   "10d-etop-high_vs_10d-etop-unsorted_lfc(2nd best)-pvalue.txt")
output_folder <- here::here("Output","Screen")
HUMAN_9606 <- read_tsv(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"),
                       col_names = FALSE) %>% set_names(c("Uniprot","Type","ID"))
Human_hsa <- inner_join( HUMAN_9606 %>% 
    subset(Type == "KEGG") %>% 
    mutate(ID = str_remove_all(ID,"hsa:")) %>% 
    dplyr::select(-Type) %>% 
        rename(KEGG = ID) %>% 
    subset(!duplicated(Uniprot)),
    HUMAN_9606 %>% 
        subset(Type == "Gene_Name") %>% 
         dplyr::select(-Type) %>% 
        subset(!duplicated(Uniprot))) %>% 
    dplyr::select(-Uniprot) %>% 
    subset(!duplicated(KEGG)) %>% 
    subset(!duplicated(ID))

Screens <- purrr::map(.x = list_of_files, ~read_tsv(here::here("Datasets","Raw",.x)) %>% 
                          mutate(lfc = lfc *(-1))) %>% 
    set_names(list_of_files)
                      
Enrichment_function <- function(Screen_input,dataset_name ){
    # Screen_input = Screens$`14d-etop_vs_14d_lfc(median)-pvalue.txt`
    # dataset_name = list_of_files[[1]]
    gene_list <- Screen_input %>% 
       left_join(Human_hsa, by = c("Gene" ="ID")) %>% 
        arrange(-lfc)
    Go_list <- gene_list %>% 
        pull(lfc,Gene)
    ego3 <- gseGO(geneList     =Go_list,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  keyType = "SYMBOL",
                  #nPerm        = 1000,
                  minGSSize    = 50,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
    if((ego3@result %>% nrow)>0){
        enrichplot::ridgeplot(ego3 %>% simplify(), showCategory = 68)+
        ggtitle(glue::glue(dataset_name," Gene BP-GSEA"))
        ggsave(here::here(output_folder, glue::glue(dataset_name," BP-GSEA.png")), height = 20, width  = 15)}
    
    ego3 <- gseGO(geneList     =Go_list,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "MF",
                  keyType = "SYMBOL",
                  #nPerm        = 1000,
                  minGSSize    = 50,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
    if((ego3@result %>% nrow)>0){
        enrichplot::ridgeplot(ego3 %>% simplify(), showCategory = 68)+
        ggtitle(glue::glue(dataset_name," Gene  MF-GSEA"))
        ggsave(here::here(output_folder, glue::glue(dataset_name," MF-GSEA.png")), height = 20, width  = 15)}
    
    KEGG_list <- gene_list%>% 
        na.omit()%>%  pull(lfc,KEGG)
    
    kk2 <- gseKEGG(geneList     = KEGG_list,
                   organism     = 'hsa',
                   minGSSize    = 10,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
    if((kk2@result %>% nrow)>0){
        
        enrichplot::ridgeplot(kk2, showCategory = 100)+
            ggtitle(glue::glue(dataset_name," Gene  KEGG "))
        ggsave(here::here(output_folder, glue::glue(dataset_name," KEGG.png")), height = 20, width  = 15)}
    
    
    mkk2 <- gseMKEGG(geneList = KEGG_list,
                     organism = 'hsa',
                     minGSSize = 10,
                     pvalueCutoff = 0.05)
    if((mkk2@result %>% nrow)>0){
        
        enrichplot::ridgeplot(mkk2, showCategory = 68)+
            ggtitle(glue::glue(dataset_name," Gene MKEGG "))
        ggsave(here::here(output_folder, glue::glue(dataset_name," MKEGG.png")), height = 20, width  = 15)}
    
    
}   
iwalk(Screens,Enrichment_function)
Screens$`14d-etop_vs_14d_lfc(median)-pvalue.txt` %>% 
    inner_join(DMSO_ruler, by = c("Gene" = "ID") ) %>% 
    ggplot(aes(x = lfc, y = Mean_DMSO, label = Gene))+
    ggrepel::geom_label_repel(data = . %>% subset(pvalue < 0.05 & !between(Mean_DMSO,-1.5,1.5))) +
    geom_point()
