#Deciding between methods 1 and 2 for DIA
#### functions used in this project ####

Load_DIA_NN_peptides <- function(report_DIA_tsv_file, Samples_df){
    read_tsv(here::here("Datasets","Raw",report_DIA_tsv_file)) %>% 
        as.data.frame() %>% 
        # subset(Lib.Q.Value<= 0.01 & Lib.PG.Q.Value <= 0.01 ) %>% 
        janitor::clean_names() %>% 
        dplyr::select(matches("run|precursor_id|precursor_normalised")) %>% 
         # distinct() %>% 
        # mutate(old_pg =protein_group) %>%  
        # subset(!str_detect(protein_group,paste0(contaminant_uniprots,collapse = "|"))) %>% 
        remove_rownames() %>%
        
        mutate(run = str_replace_all(run, c("P11690-1" = "P11591-1",
                                            "P11690-2" = "P11592-1",
                                            "P11690-3" = "P11593-1")),
               run = str_match(run,"-([:graph:]{6})-")[,2] %>% tolower()) %>%  left_join(Samples_df) %>%  
    
        pivot_wider(-run, names_from = Condition ,values_from = precursor_normalised ) %>%
    
        # group_split(protein_group) %>%
        # map_dfr(.x = .,~.x %>% discard_uniprot_isoforms("protein_group","genes" )) %>%
        # group_by(protein_group) %>%
        # mutate(Is_duplicated = n() > 1,
        #        New_Uniprot = if_else(Is_duplicated == F,protein_group,old_pg)) %>% 
        column_to_rownames("precursor_id")%>% 
        dplyr::select(any_of(Samples_df$Condition)) %>% 
        # discard_single_isoforms %>% 
        as.matrix
    
    
    
    
}
Peptide_metrics <- function(input_matrix,dataset_name){
    #this required a not normalised and not logged matrix with column names of condition in _rep format
     # dataset_name = list_files_to_analyse$DIA_report_file_1
     # input_matrix <-  Peptide_matrices$DIA_report_file_1
    
    # as.data.frame() 
    # mutate(across(everything(),~.x^2)) %>% 
    input_matrix <- input_matrix %>% as.data.frame()
    
    experimental_design_DIA <-  data.frame(
        label = colnames(input_matrix),
        condition =  str_remove_all(colnames(input_matrix),"_[:graph:]*$"),
        replicate = str_remove_all(colnames(input_matrix),"^[:graph:]*_") %>% as.numeric()
    )
    data_unique_Etop <- input_matrix %>% rownames_to_column("name") %>% 
        subset(!duplicated(name)) %>% 
        mutate(ID = name)
    
    Quant_columns <- which(colnames(data_unique_Etop) %in%colnames(input_matrix))# get LFQ column numbers
    data_se <- make_se(data_unique_Etop, Quant_columns, experimental_design_DIA)
    plot_frequency(data_se)+ggtitle(glue::glue("Peptide_overlap ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Peptide_overlap ",dataset_name,".png")))
    data_filt <- filter_missval(data_se, thr = 1)
    #data_filt2 <- filter_missval(data_se, thr = 1)
    plot_numbers(data_filt)+ggtitle(glue::glue("Peptide_numbers ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Peptide_numbers ",dataset_name,".png")))
    plot_coverage(data_filt)
    data_filt@assays@data@listData[[1]][is.nan(data_filt@assays@data@listData[[1]])] <- NA 
    png(here::here(output_folder,glue::glue("Peptide_Missingness ",dataset_name,".png")), width = 2500, height = 3800,res  =300) 
    plot_missval(data_filt)
    dev.off()
    data_norm <- normalize_vsn(data_filt)
    # data_norm@assays@data@listData[[1]] <- input_matrix
    
    DEP::meanSdPlot(data_norm)
    ggsave(here::here(output_folder,glue::glue("normalize_vsn_Peptide ",dataset_name,".png")))
    
    plot_normalization(data_se, data_norm)+ggtitle(glue::glue("Peptide_norm ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Peptide_normalisation ",dataset_name,".png")))
    
    pca_res <- prcomp(data_norm@assays@data@listData[[1]]  %>% na.omit() %>% t(), scale=TRUE)
    var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
    
    pca_res$x %>% 
        as.data.frame %>%
        rownames_to_column("Sample") %>% 
        mutate(Condition = str_remove(Sample,"_.$")) %>% 
        ggplot(aes(x=PC1,y=PC2, label = Sample, colour = Condition )) + geom_point(size=4) +
        ggrepel::geom_label_repel()+
        theme_bw(base_size=32) + 
        labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
             y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
        theme(legend.position="top") +
        ggtitle(dataset_name)+ 
        theme(plot.title = element_text(size = 20))
    ggsave(here::here(output_folder,glue::glue(dataset_name,"Peptide PCA.png")))}
Peptide_matrices <- map(.x = list_files_to_analyse[1:2],
                        ~Load_DIA_NN_peptides(.x,Samples))
walk2(Peptide_matrices,list_files_to_analyse[1:2],Peptide_metrics)
