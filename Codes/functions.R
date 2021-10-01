#### functions used in this project ####
Load_DIA_NN_Data <- function(report_DIA_tsv_file, Samples_df){
    #need to change it so that I don't end up with lines where is a single protein but the isoforms don't have ids.
    discard_single_isoforms <- function(input_matrix){
        #input_matrix = Etop_DIA
        input_matrix %>% as.data.frame() %>% 
            rownames_to_column("Uniprot") %>% 
            mutate(Uniprot_sliced = str_remove_all(Uniprot,"-[:graph:]*$")) %>% 
            group_by(Uniprot_sliced) %>% 
            mutate(Is_duplicated = n() > 1,
                   New_Uniprot = if_else(Is_duplicated == F & (!str_detect(Uniprot,";")),Uniprot_sliced,Uniprot)) %>% 
            ungroup() %>% 
            dplyr::select(-c(Is_duplicated,Uniprot_sliced,Uniprot)) %>% 
            column_to_rownames("New_Uniprot") %>% 
            as.matrix()
    }
    discard_uniprot_isoforms <- function(input_df,prot_id_col,gene_col){
        #prot_id_col = "protein_group"
        #gene_col = "genes"
        #input_df <- DIA_NN_original %>% subset(protein_group == "O60313;O60313-2")
        if(!is.na(input_df %>% pull(gene_col))){
            if(!str_detect(input_df %>% pull(gene_col),";")){
                input_df <- input_df %>% separate_rows(protein_group, sep = ";") %>%
                    left_join(HUMAN_9606_count, by = c( "protein_group"="Uniprot" )) %>%
                    arrange(n) %>% tail(1) %>% dplyr::select(-n)
            }}
        return(input_df)}
    read_tsv(here::here("Datasets","Raw",report_DIA_tsv_file)) %>% 
        as.data.frame() %>% 
        subset(Lib.Q.Value<= 0.01 & Lib.PG.Q.Value <= 0.01 ) %>% janitor::clean_names() %>% 
        dplyr::select(matches("pg_max_lfq|run|protein_group|^genes$")) %>% 
        distinct() %>% janitor::clean_names() %>% 
        mutate(old_pg =protein_group) %>%  
        subset(!str_detect(protein_group,paste0(contaminant_uniprots,collapse = "|"))) %>% 
        remove_rownames() %>%
        
        mutate(run = str_replace_all(run, c("P11690-1" = "P11591-1",
                                            "P11690-2" = "P11592-1",
                                            "P11690-3" = "P11593-1")),
               run = str_match(run,"-([:graph:]{6})-")[,2] %>% tolower()) %>%  left_join(Samples_df) %>%  
        pivot_wider(-run, names_from = Condition ,values_from = pg_max_lfq) %>% 
        group_split(protein_group) %>%
        map_dfr(.x = .,~.x %>% discard_uniprot_isoforms("protein_group","genes" )) %>%
        group_by(protein_group) %>%
        mutate(Is_duplicated = n() > 1,
               New_Uniprot = if_else(Is_duplicated == F,protein_group,old_pg)) %>% 
        column_to_rownames("New_Uniprot")%>% 
        dplyr::select(any_of(Samples_df$Condition)) %>% 
        discard_single_isoforms %>% 
        as.matrix() #%>% log2()
    
    
    
    
}

DEP_DIA <- function(input_matrix,dataset_name){
    #this required a not normalised and not logged matrix with column names of condition in _rep format
    dataset_name = "Method1_P11833_41"
    input_matrix <-  Method1_P11833_41# DIA_matrices$DIA_report_file_P11685 %>% as.data.frame()
    input_matrix <- input_matrix %>% as.data.frame()
    experimental_design_DIA <-  data.frame(
        label = colnames(input_matrix),
        condition =  str_remove_all(colnames(input_matrix),"_[:graph:]*$"),
        replicate = str_remove_all(colnames(input_matrix),"^[:graph:]*_") %>% as.numeric()
    )
    data_unique_Etop <- input_matrix %>% rownames_to_column("name") %>% 
        left_join(HUMAN_9606 %>% 
                      subset(Type == "Gene_Name") %>% 
                      dplyr::select(-Type), 
                  by = c("name" = "Uniprot")) %>%
        subset(!duplicated(name))
    
    Quant_columns <- which(colnames(data_unique_Etop) %in%colnames(input_matrix))# get LFQ column numbers
    data_se <- make_se(data_unique_Etop, Quant_columns, experimental_design_DIA)
    data_se_parsed <- make_se_parse(data_unique_Etop, Quant_columns)
    plot_frequency(data_se)+ggtitle(glue::glue("Protein_overlap ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_overlap ",dataset_name,".png")))
    data_filt <- filter_missval(data_se, thr = 1)
    #data_filt2 <- filter_missval(data_se, thr = 1)
    plot_numbers(data_filt)+ggtitle(glue::glue("Protein_numbers ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_numbers ",dataset_name,".png")))
    plot_coverage(data_filt)
    data_filt@assays@data@listData[[1]][is.nan(data_filt@assays@data@listData[[1]])] <- NA 
    png(here::here(output_folder,glue::glue("Protein_Missingness ",dataset_name,".png")), width = 2500, height = 3800,res  =300) 
    plot_missval(data_filt)
    dev.off()
    data_norm <- normalize_vsn(data_filt)
    plot_normalization(data_se, data_norm)+ggtitle(glue::glue("Protein_norm ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_normalisation ",dataset_name,".png")))
    
    # edata = data_norm@assays@data@listData[[1]] %>% na.omit()
    # batch = if_else(str_match(colnames(data_norm@assays@data@listData[[1]]),"_(.)$")[,2] %>% as.numeric() <4,2,1)
    
    # parametric adjustment
    # combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
    
    # pca_res <- prcomp(combat_edata2  %>% na.omit() %>% t(), scale=TRUE) 
    #plot_missval(data_filt)
    # Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
    
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
    ggsave(here::here(output_folder,glue::glue(dataset_name," PCA.png")))
    
    if(data_norm@assays@data@listData[[1]] %>% is.na() %>% any()){
        png(here::here(output_folder,glue::glue("Protein_Missingness_Abundance ",dataset_name,".png")), width = 2500, height = 3800,res  =300) 
        plot_detect(data_norm)
        dev.off()
        #ggsave(here::here(output_folder,glue::glue("Protein_missingness ",dataset_name,".png")))
        data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
        # to_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))>1)
        # to_not_impute <- data_norm@assays@data@listData[[1]] %>% subset(., rowSums(is.na(.))<2)
        # multiple_imputation <- impute.mi(tab = to_impute,
        #                                  conditions = experimental_design_DIA$condition %>% as.factor(),
        #                                  repbio = experimental_design_DIA$replicate %>% as.factor()) 
        # rownames(multiple_imputation) <- rownames(to_impute)
        # data_imp@assays@data@listData[[1]] <-  rbind(to_not_impute,multiple_imputation) 
    }else{
        data_imp <- data_norm
    }
    plot_imputation(data_norm, data_imp)
    ggsave(here::here(output_folder,glue::glue("Protein_imputted ",dataset_name,".png")))
    data_diff_all_contrasts <- DEP::test_diff(data_imp, type = "all")
    dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = 0)
    
    #plot_cor(dep, significant = FALSE, lower = 0, upper = 1, pal = "Reds")
    #ggsave(here::here(output_folder,glue::glue("Sample_correlation ",dataset_name,".png")))
    #plot_volcano(dep, contrast = "T0_vs_T24", label_size = 2, add_names = TRUE)
    # plot_heatmap(dep, type = "centered", kmeans = TRUE, 
    #              k = 6, col_limit = 4, show_row_names = FALSE,
    #              indicate = c("condition", "replicate"))
    # plot_single(dep, proteins = c("P05386","Q9H9B4"))
    
    data_imp@assays@data@listData[[1]] %>% 
        as.data.frame() %>% 
        rownames_to_column("ProteinGroup") %>% 
        mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
        left_join(Interesting_proteins) %>% 
        subset(!is.na(Behaviour)) %>% 
        pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        mutate(Condition= factor(Condition)) %>% 
        group_by(ProteinGroup) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>% 
        mutate(DMSO = mean(c_across(contains("dmso")), na.rm = T)) %>% 
        mutate(across(where(is.numeric), ~.x-DMSO)) %>%
        dplyr::select(-DMSO) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "ProteinGroup", "Significant")) %>% 
        ggplot(aes(x = Condition, y  = Abundance, colour = ProteinGroup,group= ProteinGroup, label = ID, alpha= Significant))+
        geom_line()+
        geom_point()+
        scale_alpha_manual(values = c(0.3,1))+
        ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DIA$label,1) & Significant == T))+
        ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DIA$label,1) & Significant == F))+
        theme(legend.position = "none") +
        facet_wrap("Behaviour")+ 
        scale_x_discrete(guide = guide_axis(n.dodge = 2))+
        ggtitle("Interesting DDR proteins Detected",
                "Significant - opaque, non-significant Transparent")
    ggsave(here::here(output_folder,glue::glue("Known_Behaviour ",dataset_name,".png")), width = 20, height = 20)
    data_norm@assays@data@listData[[1]] %>% 
        as.data.frame() %>% 
        rownames_to_column("ProteinGroup") %>% 
        mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
        pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        mutate(Condition= factor(Condition)) %>% 
        group_by(ProteinGroup) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>% 
        mutate(Samples_median = median(c_across(where(is.numeric)), na.rm = T)) %>%  
        mutate(across(where(is.numeric), ~.x-Samples_median)) %>%
        dplyr::select(-Samples_median) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "ProteinGroup", "Significant")) %>% 
        subset(Significant == T) %>%
        ggplot(aes(x = Condition, y  = Abundance, colour = ProteinGroup,group= ProteinGroup, label = ID))+
        #geom_line()+
        geom_point()+
        #scale_alpha_manual(values = c(0.3,1))+
        #ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DIA$label,1) & Significant == T))+
        #ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DIA$label,1) & Significant == F))+
        theme(legend.position = "none") +
        facet_wrap("ID")+ 
        scale_x_discrete(guide = guide_axis(n.dodge = 2))+
        ggtitle(glue::glue("Significant - Proteins",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Significant_proteins ",dataset_name,".png")), width = 20, height = 20)
    
    
    Significant_proteins <- data_imp@assays@data@listData[[1]] %>% 
        as.data.frame() %>% 
        rownames_to_column("ProteinGroup") %>% 
        #mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        #left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
        #left_join(Interesting_proteins) %>% 
        #subset(!is.na(Behaviour)) %>% 
        #pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        #mutate(Condition= factor(Condition, levels= paste(rep(c("DMSO","T0","T24"), each= 3), rep(1:3,3),sep="_"))) %>% 
        #group_by(ProteinGroup) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>% 
        #mutate(DMSO = mean(c(DMSO_1,DMSO_2,DMSO_3))) %>% mutate(across(where(is.numeric), ~.x-DMSO)) %>%
        #dplyr::select(-DMSO) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "ProteinGroup", "Significant")) %>% 
        subset(Significant == T) %>% dplyr::select(-Significant) %>% 
        pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        mutate(Condition = str_remove_all(Condition,"_.")) %>% 
        group_by(ProteinGroup,Condition) %>% 
        summarise(Mean_Abundance = mean(Abundance)) %>% 
        pivot_wider(names_from = "Condition",values_from = "Mean_Abundance") %>% 
        mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
        #mutate(duplicated = BiocGenerics::duplicated(Uniprot))
        ungroup %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>% 
        subset(!is.na(ID)) %>% 
        column_to_rownames("ID") %>%
        dplyr::select(where(is.numeric)) %>%
        mutate(across(where(is.numeric),~.x- dmso)) %>% 
        dplyr::select(-dmso) %>% 
        as.matrix() 
    paletteLength <- 50
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(Significant_proteins), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(Significant_proteins)/paletteLength, max(Significant_proteins), length.out=floor(paletteLength/2)))
    
    png(here::here(output_folder,glue::glue("Heatmap_Significant ",dataset_name,".png")), width = 2500, height = 3800,res  =300) 
    pheatmap::pheatmap(Significant_proteins,cluster_cols = F,fontsize_row = 6, clustering_distance_rows = "euclidean", 
                       main = glue::glue(dataset_name, " Significant Proteins Normalised to DMSO"),color=myColor, breaks=myBreaks)
    dev.off()
    if(names(dev.cur()) != "null device"){dev.off()}
    #clusters <- NbClust::NbClust(Significant_proteins, method = "kmeans")$Best.partition
    
    Comparisons_list <- list()
    for(i in (dep@elementMetadata %>% names() %>% str_subset("diff") )){
        i = (dep@elementMetadata %>% names() %>% str_subset("diff"))[3]
        contrast <- str_remove_all(i,"_diff")
        
        volcano_df <-  data.frame(log2_FC = dep@elementMetadata %>%  .[(glue::glue(contrast,"_diff"))] %>% unlist(),
                                  Uniprot = dep@elementMetadata$name,
                                  significant = dep@elementMetadata %>%  .[(glue::glue(contrast,"_significant"))] %>% unlist(),
                                  p.adj = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.adj"))] %>% unlist() ,
                                  p.val = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.val"))] %>% unlist()) %>% 
            mutate(Sabatini = if_else(str_detect(Uniprot,paste(Sabatini_Uniprot,collapse = "|")),T,F),
                   Single_Uniprot = Uniprot %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
            left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type), by  = c("Single_Uniprot" = "Uniprot"))
        volcano_df %>% ggplot(aes(x = log2_FC, y = -log10(p.val), label = ID, colour = Sabatini, alpha = significant))+
            geom_point()+
            ggrepel::geom_label_repel(data = . %>% subset(significant == T&  Sabatini==F ))+
            ggrepel::geom_label_repel(data = . %>% subset(significant == T&  Sabatini==T))+
            ggtitle(glue::glue("Diff Present on Chromatin", contrast),
                    subtitle = dataset_name)
        ggsave(here::here(output_folder,glue::glue("Protein_volcano_significant",dataset_name," ",contrast,".png")), width = 10, height = 15)
        
        volcano_df %>% 
            left_join(Interesting_proteins, by = c("Single_Uniprot" = "Uniprot")) %>% 
            ggplot(aes(x = log2_FC, y = -log10(p.val), label = ID, colour = Behaviour))+
            geom_point()+
            #ggrepel::geom_label_repel(data = . %>% subset(significant == T&  is.na(Behaviour)))+
            ggrepel::geom_label_repel(data = . %>% subset(!is.na(Behaviour) & p.adj<0.05), max.overlaps = 100)+
            ggtitle(glue::glue("Diff Present on Chromatin POI ", contrast),
                    subtitle = dataset_name)
        Comparisons_list[[i]] <- volcano_df
        ggsave(here::here(output_folder,glue::glue("Protein_volcano_POI",dataset_name," ",contrast,".png")), width = 10, height = 15)
        
        
        ego3 <- gseGO(geneList     = dep@elementMetadata %>% .[i] %>% unlist %>% set_names(dep@elementMetadata$name) %>% sort(decreasing = T),
                      OrgDb        = org.Hs.eg.db,
                      ont          = "MF",
                      keyType = "UNIPROT",
                      #nPerm        = 1000,
                      minGSSize    = 100,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)
        if(!is.null(ego3)){
            ridgeplot(ego3 %>% simplify(), showCategory = 68)+
                ggtitle(glue::glue(dataset_name," MF-GSEA",i))
            ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," MF-GSEA.png")), height = 20, width  = 15)}
        ego3 <- gseGO(geneList     = dep@elementMetadata %>% .[i] %>% unlist %>% set_names(dep@elementMetadata$name) %>% sort(decreasing = T),
                      OrgDb        = org.Hs.eg.db,
                      ont          = "BP",
                      keyType = "UNIPROT",
                      #nPerm        = 1000,
                      minGSSize    = 100,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)
        if(!is.null(ego3)){

            ridgeplot(ego3 %>% simplify(), showCategory = 68)+
                ggtitle(glue::glue(dataset_name,"BP-GSEA",i))
            ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," BP-GSEA.png")), height = 20, width  = 15)}
    }
    
    data_imp@assays@data@listData[[1]] %>% as.data.frame() %>%  rownames_to_column("Uniprot") %>% 
        #glimpse() %>% 
        pivot_longer(-Uniprot,names_to = "Condition", values_to = "Abundance") %>% 
        mutate(Condition = str_remove_all(Condition, "_.$")) %>% 
        group_by(Condition,Uniprot) %>% 
        summarise(Mean_Abundance = mean(Abundance, na.rm = T),
                  Protein_CV = sd(Abundance, na.rm = F)/Mean_Abundance) %>% 
        ggplot(aes(x = Mean_Abundance, y = Protein_CV))+
        geom_point()+
        facet_wrap("Condition")+
        ggtitle("Coefficient of Variation in Conditions NAs in Sd retained",
                subtitle = dataset_name)
    ggsave(here::here(output_folder,glue::glue("CV_conditions_unimputted",dataset_name," ",contrast,".png")), width = 10, height = 15)
    data_norm@assays@data@listData[[1]] %>% as.data.frame() %>%  rownames_to_column("Uniprot") %>% 
        #glimpse() %>% 
        pivot_longer(-Uniprot,names_to = "Condition", values_to = "Abundance") %>% 
        mutate(Condition = str_remove_all(Condition, "_.$")) %>% 
        group_by(Condition,Uniprot) %>% 
        summarise(Mean_Abundance = mean(Abundance, na.rm = T),
                  Protein_CV = sd(Abundance, na.rm = F)/Mean_Abundance) %>% 
        ggplot(aes(x = Mean_Abundance, y = Protein_CV))+
        geom_point()+
        facet_wrap("Condition")+
        ggtitle("Coefficient of Variation in Conditions NAs in Sd retained",
                subtitle =  dataset_name)
    ggsave(here::here(output_folder,glue::glue("CV_conditions_imputted",dataset_name," ",contrast,".png")), width = 10, height = 15)
    
    return(list(Imputted = data_imp@assays@data@listData[[1]],
                Unimputted = data_norm@assays@data@listData[[1]], 
                DEPs = set_names(Comparisons_list, dep@elementMetadata %>% names() %>% str_subset("diff"))))
}

PeCoRa_function <- function(DIA_report_file,Samples, control_condition = NULL){
    #pacman::p_load(diann, PeCorA)
    peptide_file <- read_tsv(here::here("Datasets","Raw",DIA_report_file)) 
    
    
    peptides.maxlfq <- diann_maxlfq(peptide_file[peptide_file$Lib.Q.Value <= 0.01 & peptide_file$Lib.PG.Q.Value <= 0.01,], 
                                    group.header="Stripped.Sequence",
                                    id.header = "Precursor.Id", 
                                    quantity.header = "Precursor.Normalised")
    PepCorA_setup <- peptides.maxlfq %>% as.data.frame()%>% rownames_to_column("Peptide") %>% 
        pivot_longer(-Peptide, names_to = "Condition_rep", values_to = "Normalised.Area") %>% 
        #view()
        mutate(run = str_match(Condition_rep,"-([:graph:]{6})-")[,2] %>% tolower()) %>% 
        left_join(Samples %>% mutate(run = tolower(run))) %>% #glimpse() %>% 
        separate(Condition, into = c("Condition","BioReplicate"), sep = "_") %>% 
        left_join(peptide_file %>% as.data.frame() %>% dplyr::select(Modified.Sequence,Stripped.Sequence,Protein.Group) %>% distinct(),
                  by = c("Peptide" = "Stripped.Sequence")) %>% as.data.frame() %>% mutate(Normalised.Area = replace_na(Normalised.Area,0)) %>% 
        mutate(BioReplicate = as.numeric(BioReplicate)) %>% 
        dplyr::rename(Protein=Protein.Group,
                      Peptide.Modified.Sequence = Modified.Sequence)
    area_column_name <- which(str_detect(testing_PepCorA %>% colnames(),"ormali"))
    if(is.null(control_condition)){
        control_condition <-PepCorA_setup$Condition[1]
    }
    scaled_peptides <- PeCorA_preprocessing(testing_PepCorA,
                                            area_column_name = area_column_name,
                                            threshold_to_filter = 100,
                                            control_name = control_condition)
    disagree_peptides <- PeCorA::PeCorA(scaled_peptides)
    list(disagree_peptides = disagree_peptides,
         scaled_peptides = scaled_peptides,
         disagree_proteins = disagree_peptides %>% left_join(peptide_file %>% as.data.frame() %>% dplyr::select(Genes,Protein.Group) %>% distinct(),
                                                             by = c("protein" = "Protein.Group")))
    
    
    
}
