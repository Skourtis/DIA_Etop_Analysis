#Screen
BiocManager::install("EnrichmentBrowser")
library(EnrichmentBrowser)
kegg_genes <-  EnrichmentBrowser::getGenesets("hsa",db  = "kegg")
gs <- EnrichmentBrowser::getGenesets("hsa")

library(clusterProfiler,org.Hs.eg.db)
list_new_files <- "Etop_High_low_10_cnv_negsg.mle.gene_summary_fixed_names.txt"

output_folder <- here::here("Output","Fugures")
HUMAN_9606 <- read_tsv(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"),
                       col_names = FALSE) %>% set_names(c("Uniprot","Type","ID"))
Human_hsa <- inner_join( HUMAN_9606 %>% 
    subset(Type == "KEGG") %>% 
    mutate(ID = str_remove_all(ID,"hsa:")) %>% 
    dplyr::select(-Type) %>% 
        dplyr::rename(KEGG = ID) %>% 
    subset(!duplicated(Uniprot)),
    HUMAN_9606 %>% 
        subset(Type == "Gene_Name") %>% 
         dplyr::select(-Type) %>% 
        subset(!duplicated(Uniprot))) %>% 
    dplyr::select(-Uniprot) %>% 
    subset(!duplicated(KEGG)) %>% 
    subset(!duplicated(ID))
go_KEGG <- purrr::map_dfr(.x = list(kegg_genes,gs),~
                   .x %>% unlist() %>% enframe(name = "pathway",
                                               value = "KEGG") %>% 
                   # mutate(pathway = str_remove_all(pathway,"_[:graph:]*$")) %>% 
                   separate(pathway,into = c("id","name"),sep = "_",extra = "merge") %>% 
                   left_join(Human_hsa))
fwrite(go_KEGG, here::here("Datasets","Processed","go_KEGG.csv"))
go_KEGG <- fread(here::here("Datasets","Processed","go_KEGG.csv"))
Screens_new <- purrr::map(.x = list_new_files, ~read_tsv(here::here("Datasets","Processed","Amandine_screen",.x))) %>% 
    set_names(list_new_files)
left_join(Day_10[,c(1,3)],Screens_new$Etop_High_low_10_cnv_negsg.mle.gene_summary.txt[,c(1,3)]) %>% 
    ggplot(aes(x = lfc, y = `high|beta`, label = Gene)) +
    ggrepel::geom_label_repel(data = . %>% subset( Gene %in% (Screens_new$Etop_High_low_10_cnv_negsg.mle.gene_summary.txt %>% subset(`high|fdr`<0.05) %>% pull(Gene))))+
    geom_point()+ggtitle("EtopDay10_Sort_high vs EtopDay10_unsort", subtitle = "Orginal lfc vs CNV adjusted high beta")
ggsave(here::here(output_folder,"CNV_Adjustment.png"))
Screens[["Day_10"]] <- Day_10
iwalk(.x =  Screens,~.x %>% pull(lfc) %>% hist(main = .y, breaks = 100) )                  


core_enrichment_BP <- openxlsx::read.xlsx(here::here("Datasets","Processed","BP-Multi-Omic_reduced.xlsx")) %>% 
    mutate(Genes = purrr::map(core_enrichment,str_split,pattern="/",simplify = T))

Enrichment_function <- function(Screen_input,dataset_name ){
     Screen_input =Screens_new$Etop_High_low_10_cnv_negsg.mle.gene_summary_fixed_names.txt
     dataset_name ="CNV_corrected_High_vs_unsort_EtopD10" 
     Screen_input_volc <- Screen_input %>% 
         arrange(-`high|beta`) %>% 
         mutate(Gene = as.factor(Gene),
                Rank = order(`high|beta`),
                Significant = case_when(`high|beta`<(-0.2) &`high|fdr`<0.05 ~"Low",
                                        `high|beta`>0.2& `high|fdr`<0.05 ~"High",
                                        TRUE~NA_character_))
     Screen_input_volc$Significant %>% table()
     hits_d10 <- Screen_input_volc %>% subset(!is.na(Significant)) %>% pull(Gene) 
     all_hits <- rbind(hits_d10 %>% enframe() %>% mutate(Screen = "high_g2AX"),
                   hits_d14%>% enframe() %>% mutate(Screen = "D14"))
     write_tsv(all_hits,here::here("Datasets","Processed","all_hits_CNV_cell_norm.tsv"))
     purrr::map_dfr(.x = Volcano_DFs,~.x %>% subset(  ID %in% all_hits$value) %>% dplyr::select(Uniprot,ID) ) %>% distinct %>% 
         openxlsx::write.xlsx(here::here("Datasets","Processed","Hits_Screen_in_proteome.xlsx"))
     list(high_gh2AX_screen =Screen_input_volc %>% dplyr::select(matches("Gene|gRNA|high")),
          Etoposide_survical_screen = norm_genes) %>% 
         openxlsx::write.xlsx(.,here::here("Datasets","Processed","Supp_data_1.xslx"))
     top5_hits <- Screen_input_volc %>% 
         mutate(abs_high_beta = abs(`high|beta`)) %>%
         arrange(-abs_high_beta) %>%  head(2) %>% 
         pull(`high|beta`, Gene) 
     high_low_hits <- Screen_input_volc %>% subset(!is.na(Significant)) %>% pull(Gene)
     
         ggplot(Screen_input_volc,aes(x = `high|beta`, y = Rank, label = Gene, colour = Significant))+
         geom_point()+
             ggrepel::geom_label_repel(data = . %>% subset(!is.na(Significant)),max.overlaps = 20)+
         scale_colour_manual(values = c(High = "red", Low = "blue"))+
         theme_bw()+
         ggtitle(dataset_name,
                 subtitle = "High is positive, and Low is negative selection for EtopD10_sort_high")
     ggsave(filename = here::here("Output","Figures", "CNV_corrected_High_vs_unsort_EtopD10.pdf"))
     gene_list <- Screen_input %>% 
       left_join(Human_hsa, by = c("Gene" ="ID")) %>% 
        arrange(-`high|beta`)
    Go_list <- gene_list %>% 
        pull(`high|beta`,Gene)
    
    Screen_input_volc %>%  
        mutate(Random_index= sample(Rank),
               ID = if_else(!is.na(Significant),as.character(Gene),NA_character_),
               BP = case_when(
                   Gene %in% unlist(subset(core_enrichment_KEGG, 
                                           Description == "Chemical carcinogenesis - reactive oxygen species" & comparison == "D14_BP") %>%
                                        pull(Genes))~"Chemical carcinogenesis - reactive oxygen species",
                   Gene %in% unlist(subset(core_enrichment_BP,
                                           Description == "pyrimidine deoxyribonucleotide metabolic process" & comparison == "D10_high_BP") %>%
                                        pull(Genes))~"pyrimidine deoxyribonucleotide metabolic process",
                   
                   
                   Gene %in% unlist(subset(core_enrichment_KEGG,
                                           Description == "Citrate cycle (TCA cycle)" & comparison == "D14_BP") %>%
                                        pull(Genes))~"Citrate cycle (TCA cycle)",
                   # Gene %in% unlist(subset(core_enrichment_KEGG,
                   #                         Description == "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis" & comparison == "D14_BP") %>%
                   #                      pull(Genes))~"Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
                   TRUE~NA_character_
               ),
               Significant_alpha = if_else(!is.na(Significant)|!is.na(BP),T,F)) %>% 
        ggplot(aes(x= Random_index, y= `high|beta`,label = ID,colour =BP , alpha = Significant_alpha))+
        geom_abline(slope  = 0,intercept = 0.2, alpha = 0.5,linetype="dashed")+
        geom_abline(slope  = 0,intercept = -0.2, alpha = 0.5,linetype="dashed")+
        geom_point(size = 3)+theme_bw()+
        # annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=-0.2, alpha=0.05, fill=scales::muted("red")) +
        # annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.2 , ymax=Inf, alpha=0.05, fill=scales::muted("blue")) +
        theme( panel.grid.major = element_blank(),legend.position = "none",text = element_text(size = 25),
               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        # geom_point(data = . %>% subset(is.na(Significant)), alpha = 0.1)+
         ggrepel::geom_text_repel(data = . %>% subset(!is.na(Significant)), aes(size = abs(`high|beta`)))+
        # ggrepel::geom_text_repel(size = 8)+
        
          lims(y = c(-0.51,0.51))+
    
            scale_size_continuous(range = c(7, 15)) +
        ggtitle("D10_high Difference between Etop and DMSO Annotating Significant BP terms")+
    scale_colour_manual(values = 
                            c("Chemical carcinogenesis - reactive oxygen species" = "#A95049" ,
                              "Citrate cycle (TCA cycle)" ="#363471",
                              "pyrimidine deoxyribonucleotide metabolic process" = "#1EBCD8"),
    )
              
        
    ggsave(here::here("Output","Figures","D10_rank_plot_annotated_illustrator.pdf"), height = 10, width = 20)
    D10_high_BP <- gseGO(geneList     =Go_list,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  keyType = "SYMBOL",
                  #nPerm        = 1000,
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 1,
                  #qvalueCutoff = 1,
                  verbose      = FALSE)
    # ego3 <- gseGO(geneList     =Go_list,
    #               OrgDb        = org.Hs.eg.db,
    #               ont          = "BP",
    #               keyType = "SYMBOL",
    #               #nPerm        = 1000,
    #               minGSSize    = 10,
    #               maxGSSize    = 500,
    #               pvalueCutoff = 0.05,
    #               verbose      = FALSE)
    if((D10_high_BP@result %>% nrow)>0){
        enrichplot::ridgeplot(D10_high_BP%>%  clusterProfiler::simplify() , showCategory = 68)+
        ggtitle(glue::glue(dataset_name," Gene BP-GSEA"))
        ggsave(here::here(output_folder, glue::glue(dataset_name," BP-GSEA.pdf")), height = 20, width  = 15)
        }
    
    ego3 <- gseGO(geneList     =Go_list,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "MF",
                  keyType = "SYMBOL",
                  #nPerm        = 1000,
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
    if((ego3@result %>% nrow)>0){
        enrichplot::ridgeplot(ego3 %>% clusterProfiler::simplify(), showCategory = 68)+
        ggtitle(glue::glue(dataset_name," Gene  MF-GSEA"))
        ggsave(here::here(output_folder, glue::glue(dataset_name," MF-GSEA.png")), height = 20, width  = 15)
        }
    
    KEGG_list <- gene_list%>% 
        na.omit()%>%  pull(`high|beta`,KEGG)
    
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
        ggsave(here::here(output_folder, glue::glue(dataset_name," MKEGG.png")), height = 20, width  = 15)
        }
    
    
}   
