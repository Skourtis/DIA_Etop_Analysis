#### Metabolomics Datasets #####
#Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk-16.0.2/")
# library(rJava)
# location <- getDatabase("Homo sapiens")
# pacman::p_load(knitr,cgwtools,tidyverse, MatrixCorrelation, here, ComplexHeatmap,proDA,rJava ,BridgeDbR)
# file <- "metabolites_20180508.bridge"
# download.file(
#     "https://ndownloader.figshare.com/files/10358973",
#     location
# )
# location = normalizePath(file)
# mapper <- loadDatabase(location)
# Latest features could only be obtained from this approach, make sure 'devtools' installed first
devtools::install_github("xia-lab/OptiLCMS", build = TRUE, build_vignettes = FALSE, build_manual =TRUE)
# Step 2: Install MetaboAnalystR with documentation
# devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)

library(MetaboAnalystR)
pacman::p_load("impute", "pcaMethods", "globaltest", "tidyverse",
               # "GlobalAncova", 
               "Rgraphviz", "preprocessCore",
               "genefilter",  "sva", "limma", "KEGGgraph",#"SSPA",
               "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","httr","qs")
pacman::p_load("ggplot2","parallel","reshape2","plyr",
                   "knitr","tibble","installr","fs","rmarkdown","processx","backports",
                   "bootstrap","boot","caret","dplyr","stringr","ggfortify","factoextra","MASS","FELLA",
                   "RColorBrewer","RCurl","lattice","data.table","igraph","tidyr","scales",#"MetaboAnalystR",,"DGCA"
                   "e1071","fpc","rlang","glue","digest","omu","missForest","ggpubr")
#devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)

# BiocManager::install(c("pcaMethods","limma","impute","sva","BiocParallel","genefilter","Biobase","mixOmics","statTarget", "multtest"))
# devtools::install_version("ber", version = "4.0", repos = "http://cran.us.r-project.org")
# devtools::install_version("NormalizeMets", version = "0.25", repos = "http://cran.us.r-project.org")
# devtools::install_version("metabolomics", version = "0.1.4", repos = "http://cran.us.r-project.org")
####functions####
save_pheatmap_pdf <- function(x, filename, width=1200, height=1000, res = 150) {
    pdf(filename, width = width, height = height, res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}


Savvas_anona_omu <- function (count_data, metadata, response_variable, var1, var2, 
                              interaction, log_transform, p_adjust) {
    # count_data =  Anova_test2
    #  metadata = Anova_test_Meta %>% as.data.frame() %>%
    #      mutate(across(everything(),as.factor))
    #  response_variable = "Metabolite"
    #  var1 = "Background"
    #  var2 = "Treatment"
    #  interaction = TRUE
    #  log_transform = TRUE
    #  p_adjust = "BH"
    library(plyr)
    variable1 = var1
    variable2 = var2
    rownames(count_data) <- count_data[, response_variable]
    count_data[, response_variable] <- NULL
    data_Int <- count_data[sapply(count_data, function(x) is.numeric(x))]
    data_Transpose <- as.data.frame(t(data_Int))
    data_Transpose <- as.data.frame(cbind(Sample = rownames(data_Transpose), 
                                          data_Transpose))
    data_Transpose = as.data.frame(cbind(data_Transpose, metadata))
    nums <- sapply(data_Transpose, is.numeric)
    factors <- sapply(data_Transpose, is.factor)
    data_Num <- data.frame(lapply(data_Transpose[, nums], function(x) as.numeric(x)), 
                           check.names = F, row.names = rownames(data_Transpose))
    Vect = colnames(data_Num)
    data_Ln <- log(data_Num)
    data_Ln <- as.data.frame(cbind(Sample = rownames(data_Ln), 
                                   data_Ln))
    data_Fact <- data_Transpose[, factors]
    data_Ln = merge(data_Ln, data_Fact %>% rownames_to_column("Sample") , by = "Sample")
    data_Ln <- data_Ln[, !names(data_Ln) %in% c("Sample", "Sample.1")]
    data_Num <- as.data.frame(cbind(Sample = rownames(data_Num), 
                                    data_Num))
    data_Num = merge(data_Num, data_Fact%>% rownames_to_column("Sample"), by = "Sample")
    data_Num <- data_Num[, !names(data_Num) %in% c("Sample", 
                                                   "Sample.1")]
    if (log_transform == FALSE) {
        data_mod = data_Num
    }
    else if (log_transform == TRUE) {
        data_mod = data_Ln
    }
    Mod = data_Fact
    Modifided = Mod[, !names(Mod) %in% c("Sample", "Sample.1")]
    if (missing(var2) & interaction == FALSE) {
        var1 = metadata[, var1]
        results <- llply(Vect, function(x) {
            models <- lm(data_mod[[x]] ~ var1, data = Mod)
        })
        names(results) <- Vect
        results <- lapply(results, anova)
        results <- sapply(results, cbind)
        results <- t(results)
        results <- as.data.frame(results[, "Pr(>F)"])
        results <- as.data.frame(t(results))
        colnames(results)[1] <- variable1
        colnames(results)[1] <- paste(colnames(results)[1], 
                                      "pval", sep = ".")
        results$padj = p.adjust(results[, 1], method = p_adjust)
        colnames(results)[3] <- paste(colnames(results)[3], 
                                      variable1, sep = ".")
        count_data <- cbind(rownames(count_data), data.frame(count_data, 
                                                             row.names = NULL))
        colnames(count_data)[1] <- response_variable
        results <- cbind(rownames(results), data.frame(results, 
                                                       row.names = NULL))
        colnames(results)[1] <- response_variable
        results = left_join(results, count_data, by = response_variable)
        results = results[, !(colnames(results) %in% c("V2"))]
        return(results)
    }
    else if (interaction == FALSE) {
        var1 = metadata[, var1]
        var2 = metadata[, var2]
        results <- llply(Vect, function(x) {
            models <- lm(data_mod[[x]] ~ var1 + var2, data = Mod)
        })
        names(results) <- Vect
        results <- lapply(results, anova)
        results <- sapply(results, cbind)
        results <- t(results)
        results <- as.data.frame(results[, "Pr(>F)"])
        results <- as.data.frame(t(results))
        results <- results[, 1:2]
        colnames(results)[1] <- variable1
        colnames(results)[1] <- paste(colnames(results)[1], 
                                      "pval", sep = ".")
        colnames(results)[2] <- variable2
        colnames(results)[2] <- paste(colnames(results)[2], 
                                      "pval", sep = ".")
        results$padj = p.adjust(results[, 1], method = p_adjust)
        colnames(results)[3] <- paste(colnames(results)[3], 
                                      variable1, sep = ".")
        results$padj = p.adjust(results[, 2], method = p_adjust)
        colnames(results)[4] <- paste(colnames(results)[4], 
                                      variable2, sep = ".")
        results <- cbind(rownames(results), data.frame(results, 
                                                       row.names = NULL))
        colnames(results)[1] <- response_variable
        count_data <- cbind(rownames(count_data), data.frame(count_data, 
                                                             row.names = NULL))
        colnames(count_data)[1] <- response_variable
        results = left_join(results, count_data, by = "Metabolite")
        return(results)
    }
    else if (interaction == TRUE & !missing(var2)) {
        var1 = metadata[, var1]
        var2 = metadata[, var2]
        results <- llply(Vect, function(x) {
            models <- lm(data_mod[[x]] ~ var1 + var2 + var1 * 
                             var2, data = Modifided)
        })
        names(results) <- Vect
        results <- lapply(results, anova)
        results <- sapply(results, cbind)
        results <- t(results)
        results <- as.data.frame(results[, "Pr(>F)"])
        results <- as.data.frame(t(results))
        results <- results[, 1:3]
        colnames(results)[1] <- variable1
        colnames(results)[1] <- paste(colnames(results)[1], 
                                      "pval", sep = ".")
        colnames(results)[2] <- variable2
        colnames(results)[2] <- paste(colnames(results)[2], 
                                      "pval", sep = ".")
        colnames(results)[3] <- "Interaction.pval"
        results$padj = p.adjust(results[, 1], method = p_adjust)
        colnames(results)[4] <- paste(colnames(results)[4], 
                                      variable1, sep = ".")
        results$padj = p.adjust(results[, 2], method = p_adjust)
        colnames(results)[5] <- paste(colnames(results)[5], 
                                      variable2, sep = ".")
        results$padj = p.adjust(results[, 3], method = p_adjust)
        colnames(results)[6] <- paste(colnames(results)[6], 
                                      "Interaction", sep = ".")
        results <- cbind(rownames(results), data.frame(results, 
                                                       row.names = NULL))
        colnames(results)[1] <- response_variable
        count_data <- cbind(rownames(count_data), data.frame(count_data, 
                                                             row.names = NULL))
        colnames(count_data)[1] <- response_variable
        results$padj = p.adjust(results$Interaction.pval, method = "BH")
        results = left_join(results, count_data, by = "Metabolite")
        results <- results[, !names(results) %in% c("padj")]
        return(results)
    }
}
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
mypal <- colorRampPalette( c( "blue", "grey50", "red" ) )( 50 )
output_folder <- here::here("Output","Figures")
DEP_LFQ_QC <- function(input_matrix,dataset_name){
    library(DEP)
    #this required a not normalised and not logged matrix with column names of condition in _rep format
    dataset_name = "Metabolomics"
    input_matrix <- Metabolomics_original %>%  remove_rownames() %>% 
        column_to_rownames("sample_identification") %>% 
        dplyr::select(-c(run_id,sample_description,specimen,cell_input,unit,sample_code)) %>% 
        na.omit() %>% 
        mutate(across(everything(), ~na_if(.x,"<LOD"))) %>% 
        mutate(across(everything(), ~na_if(.x,0))) %>% 
        .[,colSums(is.na(.))<(nrow(.)/2)] %>% 
        mutate(across(everything(), as.numeric)) %>% 
        t() %>% as.data.frame() #%>%
        #dplyr::select(!matches("noco|prdx1|8|24")) 
    experimental_design_DIA <-  data.frame(
        label = colnames(input_matrix),
        condition =  str_remove_all(colnames(input_matrix),"_.$"),
        replicate = str_remove_all(colnames(input_matrix),"^[:graph:]*_") %>% as.numeric()
    )
    data_unique_Etop <- input_matrix %>% rownames_to_column("name") %>% 
        left_join(Mapping %>% dplyr::select(Metabolite,KEGG), by = c("name"="Metabolite")) %>% 
        dplyr::rename(ID = KEGG)
    #### Missingness DIA SILAC LFQ###
    #### functions used in this project ####
        DEP_DIA <- function(input_matrix,dataset_name){
      input_matrix <- input_matrix %>% as.matrix
      input_matrix <- input_matrix[!(is.na(input_matrix) %>% matrixStats::rowSums2() == 24),]
      (is.na(input_matrix)) %>%  matrixStats::rowSums2() %>% set_names(rownames(input_matrix))%>% enframe(name ="Uniprot", value = "Presence") %>%
        mutate(Dataset =dataset_name,
               Presence = factor(24 - Presence,levels = 1:24))
    }
    combined <- DEP_DIA(input_matrix %>% as.data.frame() %>% dplyr::select(!matches("noco")),dataset_name)
    combined %>% ggplot(aes(x = Dataset, fill = Presence))+geom_bar(position = "stack") +
      scale_fill_manual(values =
                          c("13" = "grey99",
                            "14" = "grey95",
                            "15" = "grey90",
                            "16" ="grey80" ,
                            "18" = "grey70",
                            "19" = "grey60",
                            "20" = "grey50",
                            "21" ="grey40" ,
                            "22" = "grey30",
                            "23" = "grey20",
                            "24" = "grey10"))+
      theme_bw()+
      ggtitle("Metabolites identified and the number samples in which they were present")
 ggsave(here::here("Output","Figures","Metabolites_Missingness_qc.pdf"))   
  
    
    Quant_columns <- which(colnames(data_unique_Etop) %in%colnames(input_matrix))# get LFQ column numbers
    data_se <- make_se(data_unique_Etop, Quant_columns, experimental_design_DIA)
    #data_se_parsed <- make_se_parse(data_unique_Etop, Quant_columns)
    plot_frequency(data_se)+ggtitle(glue::glue("Protein_overlap ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_overlap ",dataset_name,".pdf")))
    data_filt <- filter_missval(data_se, thr = 1)
    #data_filt2 <- filter_missval(data_se, thr = 1)
    plot_numbers(data_filt)+ggtitle(glue::glue("Protein_numbers ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_numbers ",dataset_name,".pdf")))
    plot_coverage(data_filt)
    data_filt@assays@data@listData[[1]][is.nan(data_filt@assays@data@listData[[1]])] <- NA 
    plot_missval(data_filt)
    #data_filt@assays@data@listData[[1]] <- data_filt@assays@data@listData[[1]] %>% proDA::median_normalization()
    data_norm <- normalize_vsn(data_filt)
    #data_norm@assays@data@listData[[1]] <- data_norm@assays@data@listData[[1]] %>% #limma::normalizeMedianValues()
    #proDA::median_normalization()
    plot_normalization(data_se, data_norm)+ggtitle(glue::glue("Protein_norm ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_normalisation ",dataset_name,".pdf")))
    plot_missval(data_filt)
    # Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
    pca_res <- prcomp(data_norm@assays@data@listData[[1]] %>% as.data.frame()%>% dplyr::select(!matches("prdx|noco")) %>% na.omit() %>% t(),
                      scale=FALSE)
    var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
    
    pca_res$x %>% 
        as.data.frame %>%
        rownames_to_column("Sample") %>% 
        mutate(Condition = str_remove(Sample,"_.$")) %>% 
        ggplot(aes(x=PC1,y=PC2, label = Sample, colour = Condition )) + geom_point(size=4) +
        ggrepel::geom_label_repel()+
        theme_bw(base_size=15) + 
        labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
             y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
        theme(legend.position="right") +
        ggtitle(dataset_name)
    ggsave(here::here(output_folder,glue::glue(dataset_name," Etop PCA.pdf")),height = 15,width = 15)
    
    pca_res$rotation %>% 
      as.data.frame %>%
      rownames_to_column("Metabolite") %>% 
      left_join(DF %>% base::unclass() %>% as.data.frame()%>% dplyr::select(Metabolite, Class) ) %>% 
      #mutate(Condition = str_remove(Sample,"_.$")) %>% 
      ggplot(aes(x=PC1,y=PC2, label = Metabolite, colour = Class )) + geom_point(size=4) +
      ggrepel::geom_label_repel(data = . %>% subset(abs(PC1)>0.18 | abs(PC2)>0.18  ))+
      theme_bw(base_size=10) +
      theme(legend.position="top") +
      ggtitle(dataset_name)
    
    pca_res <- prcomp(data_norm@assays@data@listData[[1]] %>% na.omit() %>% t(),
                      scale=FALSE)
    var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
    
    pca_res$x %>% 
        as.data.frame %>%
        rownames_to_column("Sample") %>% 
        mutate(Condition = str_remove(Sample,"_.$")) %>% 
      subset(str_detect(Condition,"noco",negate = T)) %>% 
        ggplot(aes(x=PC1,y=PC2, label = Sample, colour = Condition )) + geom_point(size=4) +
        # ggrepel::geom_label_repel()+
        theme_bw(base_size=15) + 
        labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
             y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
        theme(legend.position="right") +
        ggtitle(dataset_name)
    ggsave(here::here(output_folder,glue::glue(dataset_name," PCA.pdf")),height = 15,width = 15)
    if(data_norm@assays@data@listData[[1]] %>% is.na() %>% any()){
        plot_detect(data_norm)
        ggsave(here::here(output_folder,glue::glue("Protein_missingness ",dataset_name,".pdf")))
        set.seed(1234)
        imputted <-missForest::missForest(t(data_norm@assays@data@listData[[1]]))
        data_imp <- data_norm
        data_imp@assays@data@listData[[1]] <- imputted$ximp %>% t()
        #data_imp <-     impute(data_norm, fun = "knn", rowmax = 0.8)# impute(data_norm, fun = "MinProb", q = 0.01)
            
        }else{
        data_imp <- data_norm
        }
    # write.csv(data_imp@assays@data@listData[[1]],here::here("Datasets","Processed","Metabolome_imputted.csv"))
    data_imp@assays@data@listData[[1]] <- read.csv(here::here("Datasets","Processed","Metabolome_imputted.csv"))
    plot_imputation(data_norm, data_imp)
    ggsave(here::here(output_folder,glue::glue("Protein_imputted ",dataset_name,".pdf")))
    data_diff_all_contrasts <- DEP::test_diff(data_imp, type = "manual", 
                                              test = c("wt_etop_0_vs_wt_ut", 
                                                     #  "wt_etop_8h_vs_wt_ut",
                                                       "wt_etop_24h_vs_wt_etop_0",
                                                     "wt_etop_24h_vs_wt_ut",
                                                       "wt_etop_8h_vs_wt_etop_0"  ,
                                                        "wt_etop_24h_vs_wt_etop_8h",
                                                     "prdx1_ko_ut_vs_wt_ut",
                                                     "prdx1_ko_etop_24h_vs_prdx1_ko_ut"))
    # data_diff_all_contrasts <- DEP::test_diff(data_imp, type = "manual", 
    #                                           test = c("prdx1_ko_etop_0_vs_prdx1_ko_ut", 
    #                                                    "prdx1_ko_etop_8h_vs_prdx1_ko_etop_0",
    #                                                    ))
    
    dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(0.001))
    plot_cor(dep, significant = FALSE, lower = 0, upper = 1, pal = "Reds")
    #ggsave(here::here(output_folder,glue::glue("Sample_correlation ",dataset_name,".pdf")))
    #plot_volcano(dep, contrast = "T0_vs_T24", label_size = 2, add_names = TRUE)
    # plot_heatmap(dep, type = "centered", kmeans = TRUE, 
    #              k = 6, col_limit = 4, show_row_names = FALSE,
    #              indicate = c("condition", "replicate"))
    # plot_single(dep, proteins = c("P05386","Q9H9B4"))
    
    # data_imp@assays@data@listData[[1]] %>% 
    #     as.data.frame() %>% 
    #     rownames_to_column("Metabolite") %>% 
    #     mutate(Uniprot= Metabolite %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
    #     left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
    #     left_join(Interesting_proteins) %>% 
    #     subset(!is.na(Behaviour)) %>% 
    #     pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
    #     mutate(Condition= factor(Condition, levels= paste(rep(c("DMSO","T0","T24"), each= 3), rep(1:3,3),sep="_"))) %>% 
    #     group_by(Metabolite) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>% 
    #     mutate(DMSO = mean(c(DMSO_1,DMSO_2,DMSO_3))) %>% mutate(across(where(is.numeric), ~.x-DMSO)) %>%
    #     dplyr::select(-DMSO) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
    #     left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "Metabolite", "Significant")) %>% 
    #     ggplot(aes(x = Condition, y  = Abundance, colour = Metabolite,group= Metabolite, label = ID, alpha= Significant))+
    #     geom_line()+
    #     geom_point()+
    #     scale_alpha_manual(values = c(0.3,1))+
    #     ggrepel::geom_label_repel(data = . %>% subset(Condition == "T24_3" & Significant == T))+
    #     ggrepel::geom_label_repel(data = . %>% subset(Condition == "T24_3" & Significant == F))+
    #     theme(legend.position = "none") +
    #     facet_wrap("Behaviour")+ 
    #     scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    #     ggtitle("Interesting DDR proteins Detected",
    #             "Significant - opaque, non-significant Transparent")
    # ggsave(here::here(output_folder,glue::glue("Known_Behaviour ",dataset_name,".pdf")), width = 20, height = 20)
    
    # Significant_proteins <- data_imp@assays@data@listData[[1]] %>% 
    #     as.data.frame() %>% 
    #     rownames_to_column("Metabolite") %>% 
    #     left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "Metabolite", "Significant")) %>% 
    #     subset(Significant == T) %>% dplyr::select(-Significant) %>% 
    #     #pivot_longer(-Metabolite, names_to = "Condition", values_to = "Abundance") %>% 
    #     rowwise() %>% 
    #     mutate(wt_ut= base::mean(c(wt_ut_1, wt_ut_2, wt_ut_3), na.rm = T)) %>% 
    #     #mutate(Condition = str_remove_all(Condition,"_.$")) %>% 
    #     #subset(!str_detect(Condition,"prdx|noco")) %>% 
    #     #group_by(Metabolite,Condition) %>% 
    #     #summarise(Mean_Abundance = mean(Abundance)) %>% 
    #     #pivot_wider(names_from = "Condition",values_from = "Mean_Abundance") %>% 
    #     #mutate(duplicated = BiocGenerics::duplicated(Uniprot))
    #     ungroup %>% 
    #     column_to_rownames("Metabolite") %>% 
    #     dplyr::select(where(is.numeric)) %>%
    #     mutate(across(where(is.numeric),~.x- wt_ut)) %>% 
    #     dplyr::select(-wt_ut) %>% 
    #     dplyr::select(!matches("prdx1|noco")) %>% 
    #     as.matrix() 
      Significant_proteins <- data_imp@assays@data@listData[[1]] %>% 
            as.data.frame() %>% 
            rownames_to_column("Metabolite") %>% 
            left_join(dep@elementMetadata$significant %>%
                        set_names(dep@elementMetadata$name) %>% enframe(name = "Metabolite", "Significant")) %>% 
            subset(Significant == T) %>% dplyr::select(-Significant) %>% 
            pivot_longer(-Metabolite, names_to = "Condition", values_to = "Abundance") %>% 
            mutate(Condition = str_remove_all(Condition,"_.$")) %>% 
            subset(!str_detect(Condition,"prdx|noco")) %>% 
            group_by(Metabolite,Condition) %>% 
            dplyr::summarise(Mean_Abundance = mean(Abundance)) %>% 
            pivot_wider(names_from = "Condition",values_from = "Mean_Abundance") %>% 
            mutate(duplicated = BiocGenerics::duplicated(Metabolite)) %>% 
        ungroup %>% 
            column_to_rownames("Metabolite") %>% 
            dplyr::select(where(is.numeric)) %>%
            mutate(across(where(is.numeric),~.x- wt_ut)) %>% 
            dplyr::select(-wt_ut) %>% 
            as.matrix() 
    paletteLength <- 50
    myColor <- colorRampPalette(c("#1B345D", "white", "#6E1927"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(Significant_proteins), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(Significant_proteins)/paletteLength, max(Significant_proteins), length.out=floor(paletteLength/2)))
    myBreaks <- seq(-2.5,2.5, length.out= 50)
    annotation_row= data.frame(name = rownames(Significant_proteins)) %>% 
        left_join(dep@elementMetadata %>% as.data.frame() %>% dplyr::select(matches("name|_significant"))) %>% 
        column_to_rownames("name") %>% 
        mutate(across(everything(),~if_else(. ==T,"Significant","Not_Significant"))) %>% 
        set_names(.,str_remove_all(names(.),"_significant"))%>% 
        .[rev(c(3,7,5))] %>%
      rownames_to_column("Metabolite") 
    pathway <- data.table(
      Metabolite =  annotation_row$Metabolite,
               Pathway = c("Nucleotide metabolism","Citrate cycle (TCA cycle)","Arginine biosynthesis","Nucleotide metabolism",
                           "Nucleotide metabolism","Pentose and glucuronate interconversion","Pentose phosphate pathway","Nucleotide metabolism",
                           "Nucleotide metabolism","Nucleotide metabolism","Pyruvate metabolism","Tyrosine metabolism",
                           "Arginine biosynthesis","Nucleotide metabolism","Glycerophospholipid metabolism","Riboflavin metabolism",
                           "Pentose phosphate pathway","Nucleotide metabolism","Pentose phosphate pathway"
                           ))
    annotation_row <- left_join(annotation_row,pathway)  %>% 
      add_count(Pathway) %>% arrange(-n) %>% dplyr::select(-n)%>%  column_to_rownames("Metabolite") 
    my_colour = rep(list(c(Significant = "#011480", Not_Significant = "grey90")),length(colnames(annotation_row))) %>% 
        set_names(colnames(annotation_row))
    set.seed(124)
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    my_colour$Pathway <-sample(color, length(annotation_row$Pathway %>% unique())) %>% set_names(annotation_row$Pathway %>% unique())
    pdf(here::here(output_folder,glue::glue("Heatmap_Significant ",dataset_name,".pdf")),width = 11,height =12)
    Significant_Etop <- pheatmap::pheatmap(Significant_proteins[rownames(annotation_row),c(1,3,2)],cluster_cols = F,
                                          cluster_rows = F, fontsize =14,fontsize_col = 12,fontsize_row = 13,
                       annotation_row = annotation_row,annotation_colors = my_colour,
                       main = glue::glue(dataset_name, " Significant Normalised to Untreated"),color=myColor, breaks=myBreaks)
    dev.off()
    Significant_proteins <- data_imp@assays@data@listData[[1]] %>%
      as.data.frame() %>%
      rownames_to_column("Metabolite") %>%
      left_join(dep@elementMetadata$significant %>%
                  set_names(dep@elementMetadata$name) %>% enframe(name = "Metabolite", "Significant")) %>%
      # subset(Significant == T) %>%
      dplyr::select(-Significant) %>%
      pivot_longer(-Metabolite, names_to = "Condition", values_to = "Abundance") %>%
      mutate(Condition = str_remove_all(Condition,"_.$")) %>%
      subset(!str_detect(Condition,"prdx|noco")) %>%
      group_by(Metabolite,Condition) %>%
      dplyr::summarise(Mean_Abundance = mean(Abundance)) %>%
      pivot_wider(names_from = "Condition",values_from = "Mean_Abundance") %>%
      mutate(duplicated = BiocGenerics::duplicated(Metabolite)) %>%
      ungroup %>%
      column_to_rownames("Metabolite") %>%
      dplyr::select(where(is.numeric)) %>%
      mutate(across(where(is.numeric),~.x- wt_ut)) %>%
      dplyr::select(-wt_ut) %>%
      as.matrix()
    paletteLength <- 50
    myColor <- colorRampPalette(c("#1B345D", "white", "#6E1927"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(Significant_proteins), 0, length.out=ceiling(paletteLength/2) + 1),
                  seq(max(Significant_proteins)/paletteLength, max(Significant_proteins), length.out=floor(paletteLength/2)))
    myBreaks <- seq(-2.5,2.5, length.out= 50)
    pdf(here::here(output_folder,glue::glue("Heatmap_All_metabolites ",dataset_name,".pdf")),width = 11,height =12)
    Significant_Etop <- pheatmap::pheatmap(Significant_proteins[,c(1,3,2)],cluster_cols = F,
                                           cluster_rows = F, fontsize =14,fontsize_col = 12,fontsize_row = 13,
                                           # annotation_row = annotation_row,annotation_colors = my_colour,
                                           main = glue::glue(dataset_name, " All_metabo Normalised to Untreated"),color=myColor, breaks=myBreaks)


    dev.off()
    
    Significant_proteins <- data_imp@assays@data@listData[[1]] %>%
      as.data.frame() %>%
      rownames_to_column("Metabolite") %>%
      left_join(dep@elementMetadata$significant %>%
                  set_names(dep@elementMetadata$name) %>% enframe(name = "Metabolite", "Significant")) %>%
      # subset(Significant == T) %>%
      dplyr::select(-Significant) %>%
      pivot_longer(-Metabolite, names_to = "Condition", values_to = "Abundance") %>%
      mutate(Condition = str_remove_all(Condition,"_.$")) %>%
      subset(!str_detect(Condition,"wt|noco")) %>%
      group_by(Metabolite,Condition) %>%
      dplyr::summarise(Mean_Abundance = mean(Abundance)) %>%
      pivot_wider(names_from = "Condition",values_from = "Mean_Abundance") %>%
      mutate(duplicated = BiocGenerics::duplicated(Metabolite)) %>%
      ungroup %>%
      column_to_rownames("Metabolite") %>%
      dplyr::select(where(is.numeric)) %>%
      mutate(across(where(is.numeric),~.x- prdx1_ko_ut)) %>%
      dplyr::select(-prdx1_ko_ut) %>%
      as.matrix()
    paletteLength <- 50
    myColor <- colorRampPalette(c("#1B345D", "white", "#6E1927"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(Significant_proteins), 0, length.out=ceiling(paletteLength/2) + 1),
                  seq(max(Significant_proteins)/paletteLength, max(Significant_proteins), length.out=floor(paletteLength/2)))
    myBreaks <- seq(-2.5,2.5, length.out= 50)
    pdf(here::here(output_folder,glue::glue("Heatmap_All_metabolites_PRDX1 ",dataset_name,".pdf")),width = 11,height =12)
    Significant_Etop <- pheatmap::pheatmap(Significant_proteins[,c(1,3,2)],cluster_cols = F,
                                           cluster_rows = F, fontsize =14,fontsize_col = 12,fontsize_row = 13,
                                           # annotation_row = annotation_row,annotation_colors = my_colour,
                                           main = glue::glue(dataset_name, " All_metabo Normalised to Untreated"),color=myColor, breaks=myBreaks)
    
    
    dev.off()
    
    # save_pheatmap_pdf(Significant_Etop, here::here(output_folder, "my_heatmap.pdf"),width=1600, height=800, )
    #clusters <- NbClust::NbClust(Significant_proteins, method = "kmeans")$Best.partition
    
     Comparisons_list <- list()
    for(i in (dep@elementMetadata %>% names() %>% str_subset("diff") )){
          # i = (dep@elementMetadata %>% names() %>% str_subset("diff")) %>% .[1]
        contrast <- str_remove_all(i,"_diff")
        
        volcano_df <-  data.frame(log2_FC = dep@elementMetadata %>%  .[(glue::glue(contrast,"_diff"))] %>% unlist(),
                                  Metabolite = dep@elementMetadata$name,
                                  significant = dep@elementMetadata %>%  .[(glue::glue(contrast,"_significant"))] %>% unlist(),
                                  p.adj = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.adj"))] %>% unlist() ,
                                  p.val = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.val"))] %>% unlist()) %>% 
            left_join(Type)
           
        # volcano_df %>% 
        #     ggplot(aes(x = log2_FC, y = -log10(p.val), label = Metabolite,  alpha = significant, colour = Class))+
        #     geom_point()+
        #     ggrepel::geom_label_repel(data = . %>% subset(significant == T), max.overlaps = 200)+
        #     #ggrepel::geom_label_repel(data = . %>% subset(significant == T&  Sabatini==T))+
        #     ggtitle(glue::glue("Diff abundance ", contrast))+
        #   theme_bw()+
        #     xlim(c(-4,4))+
        #     ylim(c(0,20))
        Comparisons_list[[i]] <- volcano_df
         # ggsave(here::here(output_folder,glue::glue("Metabo_volcano_significant",contrast,".pdf")), width = 8, height = 8)
    }
     names(Comparisons_list) <- names(Comparisons_list) %>% str_remove_all("prdx1_")
    openxlsx::write.xlsx(Comparisons_list,here::here("Datasets","Processed","DEPS_metabolomics.xlsx"),overwrite = T)
}
###### Data loading#####


Metabolomics_original <- openxlsx::read.xlsx(here::here("Datasets","Raw", "Report_MCC077_METAB02_20210728.xlsx"),startRow = 3 ) %>% 
    janitor::clean_names() %>% 
    mutate(sample_identification = janitor::make_clean_names(sample_identification)) %>% 
    subset(!is.na(sample_code))
Conditions <- Metabolomics_original %>% dplyr::select(sample_description,sample_identification) %>% 
    na.omit() %>% mutate(sample_identification = sample_identification %>%  str_remove("_.$")) %>% distinct()
# testing_normalisation <- Metabolomics_original %>%  remove_rownames() %>% 
#     column_to_rownames("sample_identification") %>% 
#     dplyr::select(-c(run_id,sample_description,specimen,unit,sample_code)) %>% 
#     na.omit() %>% 
#     mutate(across(everything(), ~na_if(.x,"<LOD"))) %>% 
#     mutate(across(everything(), ~na_if(.x,0))) %>% 
#     .[,colSums(is.na(.))<(nrow(.)/2)] %>% 
#     mutate(across(everything(), as.numeric),
#            across(-cell_input, log2)) %>% 
#     mutate(cell_input = cell_input/100000,
#            across(everything(),~.x/cell_input )) %>% 
#     dplyr::select(-cell_input) #%>% as.matrix%>% t() %>% proDA::median_normalization()
# testing_normalisation %>% t() %>%  boxplot()
Metabolomics <-  
    Metabolomics_original %>%  remove_rownames() %>% 
    column_to_rownames("sample_identification") %>% 
    dplyr::select(-c(run_id,sample_description,specimen,cell_input,unit,sample_code)) %>% 
    na.omit() %>% 
    mutate(across(everything(), ~na_if(.x,"<LOD"))) %>% 
    mutate(across(everything(), ~na_if(.x,0))) %>% 
    .[,colSums(is.na(.))<(nrow(.)/2)] %>% 
    mutate(across(everything(), as.numeric)) %>% 
    t()%>% proDA::median_normalization()# %>% log2()
Metabolomics <-  read.csv(here::here("Datasets","Processed","Metabolome_imputted.csv")) %>% 
  as.data.frame() %>% column_to_rownames("X")
Mapping <- data.frame(Metabolite= rownames(Metabolomics),
                          Query = rownames(Metabolomics)%>% 
               str_replace_all("_"," ") %>%
               str_replace_all(" 6 ", "-6-") %>% 
               str_replace_all(" 1 ", "-1-")%>%
               str_replace_all(" 5 ", "-5-")%>% 
               str_replace_all("^n ", "n-") %>%
               str_replace_all("^o ", "o-")) %>% 
    left_join(readr::read_csv(here::here("Datasets","Raw","Metabolite_MetaboAnalyst_name_mapping.csv")) %>% 
    dplyr::select(Query,KEGG) %>% na.omit() %>% 
    rbind(data.frame(Query = c("pyridoxal hydrochloride","hydroxy glutamic acid","pentahydroxyhexanoic acid","oxamic acid","hexose-6-phosphate",
                     "phospho serine","gamma glu cys","hexose-1-phosphate","dihydroxyisovalerate","deoxy methylthio adenosine",
                     "ketobutyrate","gluthathione oxidized","n-carbamoyl aspartic acid","n-carbamyl glutamic acid","hydroxyglutaric acid",
                     "uridine diphosphohexose","alpha ketoglutaric acid","phosphoglyceric acid","fructose 1-6-biphosphate","pyridoxal-5-phosphate",
                     "citric acid isocitric acid","phenolred","ketovaleric acid","myristoylcarnitine"),
           KEGG = c("C00250","C03079","C00257","C01444","C02965",
                    "C01005","C00669","C01002","C04039","C00170",
                    "C00109","C00127","C00438","C05829","C02630",
                    "C00029","C00026","C00597","C00354","C00018",
                    "C00311","C12600",'C06255',"C00487"))))
testing <- Metabolomics%>% #log2() %>% 
   # impute::impute.knn() %>% .[["data"]] %>% 
    as.data.frame() %>% 
    rownames_to_column("Metabolite") %>%    rowwise() %>% 
    mutate(Ctrl= base::mean(c(wt_ut_1, wt_ut_2, wt_ut_3), na.rm = T)) %>% 
   ungroup() %>% 
    mutate(across(where(is.numeric),~.x-Ctrl))%>% dplyr::select(-Ctrl) %>% 
    pivot_longer(-Metabolite, names_to = "Sample",values_to = "Abundance") %>% 
    mutate(Sample = Sample %>% str_replace_all("prdx1_ko","prdx1ko") %>% 
               str_replace_all("ut_","ut_wt_" )) %>% 
    
    separate(Sample, c("Background","Treatment","Time","Replicate")) %>% 
    mutate(Time = if_else(Treatment == "ut","-3",Time)) %>% 
    mutate(Time = factor(Time, levels= c("-3","0","8h","24h"))) %>% 
    subset(Treatment != 'noco')
# testing_MetaboAnalyst <- Metabolomics%>% 
#     impute::impute.knn() %>% .[["data"]] %>% as.data.frame() %>% 
#     rownames_to_column("Metabolite") %>%    #rowwise() %>% 
#     #mutate(Ctrl= base::mean(c(wt_ut_1, wt_ut_2, wt_ut_3), na.rm = T)) %>% 
#     #ungroup() %>% 
#     #mutate(across(where(is.numeric),~.x-Ctrl))%>% dplyr::select(-Ctrl) %>% 
#     pivot_longer(-Metabolite, names_to = "Sample",values_to = "Abundance") %>% 
#     mutate(Sample = Sample %>% str_replace_all("prdx1_ko","prdx1ko") %>% 
#                str_replace_all("ut_","ut_-3_" ) ,
#            Phenotype =Sample %>% str_match("^([:graph:]*?)_") %>% .[,2],
#            Time =Sample %>% str_match("_[:graph:]*?_([:graph:]*?)_.$") %>% .[,2] %>% str_remove_all("h") %>% as.numeric) %>% 
#     pivot_wider(names_from = "Metabolite",values_from = "Abundance") #%>% t()
#     separate(Sample, c("Background","Treatment","Time","Replicate")) %>% 
#     mutate(Time = factor(Time)) %>% 
#     subset(Treatment != 'noco')
# testing_noreva2020 <- Metabolomics%>% as.data.frame() %>% 
#       rownames_to_column("Metabolite") %>%    
#     pivot_longer(-Metabolite, names_to = "Sample",values_to = "Abundance") %>% 
#     mutate(Sample = Sample %>% str_replace_all("prdx1_ko","prdx1ko") %>% 
#                str_replace_all("ut_","ut_-3_" ) ,
#            Phenotype =Sample %>% str_match("^([:graph:]*?)_") %>% .[,2],
#            Time =Sample %>% str_match("_[:graph:]*?_([:graph:]*?)_.$") %>% .[,2] %>% str_remove_all("h") %>% as.numeric) %>% 
#     pivot_wider(names_from = "Metabolite",values_from = "Abundance") %>% #t()
#     rename(label= Phenotype, time = Time) %>% 
#     subset(!str_detect(Sample,"noco")) %>% 
#     mutate(Replicate = str_match(Sample,"(.)$") %>% .[,2],
#            time := glue::glue("T{time+3}")) %>% 
#     rename(Phenotype=label,
#            Time = time) %>% 
#     mutate(Time = case_when(
#         Time == "T0" ~"T0",
#         Time == "T3"~"T1",
#         Time == "T11"~"T2",
#         Time == "T27"~"T3"),
#         Sample = glue::glue("S{Replicate}{Time}")) %>% 
#     dplyr::select(-c(Replicate)) 
# write.csv(testing_noreva2020 %>% subset(Phenotype =="wt"),here::here("Project_Output","testing_noreva2020.csv"),row.names = F)


Metabolomics_count <- Metabolomics%>% 
    #impute::impute.knn() %>% .[["data"]] %>% 
    as.data.frame()%>% rownames_to_column("Metabolite") %>% left_join(Mapping) %>% 
    dplyr::select(-Query) 
Meta_data <- Metabolomics_original %>% 
    dplyr::select(sample_code,sample_identification) %>%
    mutate(Sample = sample_identification,
           sample_identification = str_remove_all(sample_identification,"_.$") %>% 
               str_replace_all("prdx1_ko","prdx1ko") %>% 
               str_remove_all("h")) %>% 
    separate(sample_identification,c("Background","Treatment","Time")) %>% 
    mutate(Time = replace_na(Time, 0) %>% as.numeric) %>% 
    na.omit() %>% 
    #subset(Time == 0) %>% 
    #subset(Background =="wt") %>% 
    subset(Treatment != "noco") %>% 
    dplyr::select(-Time) %>% 
    rowwise() %>% 
    mutate(Grouped = paste0(Background,Treatment)) %>% ungroup() %>% 
    mutate(Background = as.factor(Background),
           Treatment = as.factor(Treatment)
           )
DF <- omu::assign_hierarchy(count_data = Metabolomics_count%>% dplyr::select(!matches("noco|prdx1")) , keep_unknowns = TRUE, identifier = "KEGG")
Type <- DF[,c("Metabolite","KEGG","Class")] %>% base::unclass() %>% as.data.frame() %>% 
  mutate(Class = as.character(Class), Class = if_else(KEGG == "C00081","Nucleic acids",Class))
#DF <- assign_hierarchy(count_data = c57_nos2KO_mouse_countDF, keep_unknowns = TRUE, identifier = "KEGG")
excluded <- Metabolomics_original %>% 
    subset(!(sample_code %in% Meta_data$sample_code)) %>% pull(sample_identification)
Anova_test <- Metabolomics_count %>% 
    dplyr::select(-excluded) %>% 
    mutate(across(where(is.numeric),exp))
Anova_test_Meta <-Meta_data %>% 
    dplyr::select(Sample, Background, Treatment, Grouped)
write.csv(Anova_test, here::here("Datasets","Processed","Anova_omu_test.csv"), row.names = F)
Anova_test2 <- omu::read_metabo(filepath = here::here("Datasets","Processed","Anova_omu_test.csv")) 
# %>% 
#     .[,c("Metabolite", "KEGG"    , "wt_ut_1" ,          "wt_ut_2"  ,         "wt_ut_3",           "wt_etop_0_1" ,
#          "wt_etop_0_2" ,      "wt_etop_0_3"    ,"prdx1_ko_ut_1"   ,  "prdx1_ko_ut_2"    , "prdx1_ko_ut_3"   ,
#          "prdx1_ko_etop_0_1" ,"prdx1_ko_etop_0_2", "prdx1_ko_etop_0_3")]
savvas_omu <- function (count_data, metadata, numerator, denominator, response_variable, 
                        Factor, log_transform, p_adjust, test_type) 
{
    count_data = DF_one_factor
    metadata = Meta_data %>% as.data.frame() %>%subset(Background == "wt" & str_detect(Sample, "wt_etop_0|wt_ut")) %>% mutate(Grouped = as.factor(Grouped))
    numerator = "etop"
    denominator = "ut"
    response_variable = "Metabolite"
    Factor = "Treatment"
    log_transform = TRUE
    p_adjust = "fdr"
    test_type = "welch"
    
    rownames(count_data) <- count_data[, response_variable]
    count_data[, response_variable] <- NULL
    data_Int <- count_data[sapply(count_data, function(x) is.numeric(x))]
    data_Transpose <- as.data.frame(t(data_Int))
    data_Transpose <- as.data.frame(cbind(Sample = rownames(data_Transpose), 
                                          data_Transpose))
    Factor = metadata[, Factor]
    data_Transpose$Factor = Factor[match(metadata$Sample, data_Transpose$Sample)]
    data_Subset <- filter(data_Transpose, Factor == numerator | 
                              Factor == denominator)
    rownames(data_Subset) <- data_Subset[, "Sample"]
    data_Subset[, "Sample"] <- NULL
    data_Numeric <- data_Subset[sapply(data_Subset, function(x) is.numeric(x))]
    data_Numeric <- data.frame(lapply(data_Numeric, function(x) as.numeric(as.character(x))), 
                               check.names = F, row.names = rownames(data_Numeric))
    data_Log <- as.data.frame(log(data_Numeric))
    data_Log$Factor = data_Subset$Factor
    cols_to_test <- data_Log[sapply(data_Log, function(x) is.numeric(x))]
    Vect = colnames(cols_to_test)
    data_Numeric$Factor <- data_Subset$Factor
    if (log_transform == FALSE) {
        data_mod = data_Numeric  }
    else if (log_transform == TRUE) {
        data_mod = data_Log
    }
    model = data_mod[, "Factor"]
    if (test_type == "students") {
        Run_Tests <- function(data_Log, Vect, model) {
            results <- ldply(Vect, function(Metabolite) {
                t_val = t.test(data_mod[[Metabolite]] ~ model)$statistic
                p_val = t.test(data_mod[[Metabolite]] ~ model)$p.value
                return(data.frame(Metabolite = Metabolite, t_value = t_val, 
                                  pval = p_val))
            })
        }
    }
    else if (test_type == "mwu") {
        Run_Tests <- function(data_Log, Vect, model) {
            results <- ldply(Vect, function(Metabolite) {
                t_val = wilcox.test(data_mod[[Metabolite]] ~ 
                                        model)$statistic
                p_val = wilcox.test(data_mod[[Metabolite]] ~ 
                                        model)$p.value
                return(data.frame(Metabolite = Metabolite, t_value = t_val, 
                                  pval = p_val))
            })
        }
    }
    else if (test_type == "welch") {
        Run_Tests <- function(data_Log, Vect, model) {
             # data_Log = data_Log; Vect = Vect; model = model
            results <- plyr::ldply(Vect, function(Metabolite) {
                 # Metabolite = Vect[1]
                t_val = t.test(data_mod[[Metabolite]] ~ model, 
                               var.equal = FALSE)$statistic
                p_val = t.test(data_mod[[Metabolite]] ~ model, 
                               var.equal = FALSE)$p.value
                return(data.frame(Metabolite = Metabolite, t_value = t_val, 
                                  pval = p_val))
            })
        }
    }
    results <- Run_Tests(data_Log = data_Log, Vect = Vect, model = model)
    data_Log$Factor <- factor(data_Log$Factor, levels = c(numerator, 
                                                          denominator))
    data_Numeric$Factor = data_Log$Factor
    Means <- data_Numeric %>% group_by(Factor) %>% summarise_all(funs(mean))
    Means_T = as.data.frame(t(Means))
    colnames(Means_T) <- as.character(unlist(Means_T[1, ]))
    Means_T = Means_T[-1, ]
    colnames(Means_T)[1] <- numerator
    colnames(Means_T)[1] <- paste(colnames(Means_T)[1], "mean", 
                                  sep = ".")
    colnames(Means_T)[2] <- denominator
    colnames(Means_T)[2] <- paste(colnames(Means_T)[2], "mean", 
                                  sep = ".")
    Means_T[, 3] <- as.numeric(as.character(Means_T[, 1]))/as.numeric(as.character(Means_T[, 
                                                                                           2]))
    colnames(Means_T)[3] <- "Fold_Change"
    Means_T[, 4] = log2(as.numeric(as.character(Means_T[, 3])))
    colnames(Means_T)[4] = "log2FoldChange"
    Means_T <- cbind(rownames(Means_T), data.frame(Means_T, 
                                                   row.names = NULL))
    colnames(Means_T)[1] <- response_variable
    st.err <- function(x) sd(x)/sqrt(length(x))
    stdev <- data_Log %>% group_by(Factor) %>% summarise_all(funs(sd))
    stdev <- as.data.frame(stdev)
    rownames(stdev) <- stdev[, "Factor"]
    stdev[, "Factor"] <- NULL
    stdev_t = as.data.frame(t(stdev))
    stdev_t <- cbind(rownames(stdev_t), data.frame(stdev_t, 
                                                   row.names = NULL))
    colnames(stdev_t)[1] <- response_variable
    colnames(stdev_t)[2] <- numerator
    colnames(stdev_t)[2] <- paste(colnames(stdev_t)[2], "stdev", 
                                  sep = ".")
    colnames(stdev_t)[3] <- denominator
    colnames(stdev_t)[3] <- paste(colnames(stdev_t)[3], "stdev", 
                                  sep = ".")
    SE <- data_Log %>% group_by(Factor) %>% summarise_all(funs(st.err))
    SE <- as.data.frame(SE)
    rownames(SE) <- SE[, "Factor"]
    SE[, "Factor"] <- NULL
    SE_t = as.data.frame(t(SE))
    SE_t <- cbind(rownames(SE_t), data.frame(SE_t, row.names = NULL))
    colnames(SE_t)[1] <- response_variable
    colnames(SE_t)[2] <- numerator
    colnames(SE_t)[2] <- paste(colnames(SE_t)[2], "std.err", 
                               sep = ".")
    colnames(SE_t)[3] <- denominator
    colnames(SE_t)[3] <- paste(colnames(SE_t)[3], "std.err", 
                               sep = ".")
    colnames(results)[1] <- response_variable
    results$padj = p.adjust(results$pval, method = p_adjust)
    count_data = cbind(rownames(count_data), data.frame(count_data, 
                                                        row.names = NULL))
    colnames(count_data)[1] <- "Metabolite"
    results = merge(results, count_data, by = response_variable, 
                    all = TRUE)
    results = merge(Means_T, results, by = response_variable, 
                    all.y = TRUE)
    results = merge(SE_t, results, by = response_variable, all.y = TRUE)
    results = merge(stdev_t, results, by = response_variable, 
                    all.y = TRUE)
    class(results) = append(class(results), "cpd")
    return(results)
}

# DF_one_factor <- omu::assign_hierarchy(count_data = Metabolomics_count%>% dplyr::select(!matches("noco|prdx1|8|24")) %>% 
#                                       mutate(across(where(is.numeric),exp)), keep_unknowns = TRUE, identifier = "KEGG")
# 
# DF_stats <- omu::omu_summary(count_data = DF_one_factor, metadata = Meta_data %>% as.data.frame() %>% 
#                             subset(Background == "wt" & str_detect(Sample, "wt_etop_0|wt_ut")), 
#                         numerator = "etop", denominator = "ut", response_variable = "Metabolite",
#                         Factor = "Treatment", log_transform = TRUE, p_adjust = "fdr", test_type = "welch")


# DF_stats_counts <- count_fold_changes(count_data = DF_stats, "Class", column = "Class", sig_threshold = 0.05, keep_unknowns = T)
# library(ggplot2)
# plot_bar(fc_data = DF_stats_counts, fill = c("dodgerblue2", "firebrick2"), color = c("black", "black"), size = c(1,1)) + labs(x = "Class") + theme(panel.grid = element_blank())
# DF_stats %>% ggplot(aes(x = log2FoldChange, y = -log10(pval),colour = Class, label = Metabolite))+
#     geom_point()+
#     ggrepel::geom_label_repel(data =. %>% subset(padj<0.05))
# openxlsx::write.xlsx(DF_stats %>% base::unclass() %>% as.data.frame() %>% dplyr::select(KEGG,log2FoldChange) %>% 
#                          mutate(Colour = map2color(log2FoldChange,mypal),
#                                 Width = "W22"), here::here("Project_Output","testingipath3.xlsx"))
# 
# DF_stats %>% ggplot(aes(x= log2FoldChange, y = -log10(pval), colour = Class, label = Metabolite))+
#     geom_point()+
#     ggrepel::geom_text_repel(data = . %>% subset(padj<0.05))
# DF_anova <- Savvas_anona_omu(count_data =  Anova_test2 %>% dplyr::select(-contains("prdx")), 
#                              metadata = Anova_test_Meta %>% as.data.frame() %>% 
#                                  mutate(across(everything(),as.factor)) %>% 
#                                  subset(str_detect(Background,"prdx")),
#                              response_variable = "Metabolite",
#                              var1 = "Treatment",
#                              #var2 = "Treatment",
#                              interaction = FALSE,
#                              log_transform = FALSE,
#                              p_adjust = "BH") # %>%  
library(plyr)
DF_anova <- Savvas_anona_omu(count_data =  Anova_test2, 
                             metadata = Anova_test_Meta %>% as.data.frame() %>% 
                                 mutate(across(everything(),as.factor)),
                             response_variable = "Metabolite",
                             var1 = "Background",
                             var2 = "Treatment",
                             interaction = TRUE,
                             log_transform = TRUE,
                             p_adjust = "BH")  %>%  
    mutate(Significant_type= case_when(
        padj.Background<0.05 & padj.Treatment<0.05 ~ "Significant Background and Treatment",
        padj.Background<0.05 ~ "Significant Background",
        padj.Treatment<0.05 ~ "Significant Treatment",
        padj.Interaction<0.05 ~ "Significant Interaction",
        TRUE~"Not Significant")) %>% 
       mutate(Significant = if_else(Significant_type != "Not Significant",T,F))
 DF_anova %>%
   left_join(Type) %>% 
   mutate(Class = if_else(Class %in% c("Nucleic acids","Vitamins and Cofactors"),Class,NA_character_),
          Metabolite = str_replace_all(Metabolite,"_","-")) %>% 
   ggplot(aes(x = -log10(Background.pval), y = -log10(Treatment.pval),colour = Class,  label = Metabolite))+
    geom_point(aes(alpha =Significant), stroke = 0.1,size = 3)+
    ggrepel::geom_text_repel(data = . %>% subset(Significant ==T & Significant_type == "Significant Treatment" & (-log10(Background.pval))<0.8),fontface = 'bold',  aes(colour = Class),size = 8)+
   ggrepel::geom_text_repel(data = . %>% subset(Significant ==T & Significant_type == "Significant Background and Treatment" & (-log10(Background.pval))>6),fontface = 'bold', aes(colour = Class),size = 8)+
   ggrepel::geom_text_repel(data = . %>% subset(Significant ==T & Significant_type == "Significant Background" & (-log10(Treatment.pval))<1.5 & (-log10(Background.pval))>3),fontface = 'bold', aes(colour = Class),size = 8)+
   annotate("rect", xmin=0 , xmax=-log10(0.01), ymin=-log10(0.01), ymax=Inf, alpha=0.05, fill=scales::muted("red")) +
   annotate("rect", xmin=-log10(0.01) , xmax=Inf, ymin=0, ymax=-log10(0.01), alpha=0.05, fill=scales::muted("blue")) +
   # annotate("rect", xmin=-log10(0.01) , xmax=Inf, ymin=-log10(0.01), ymax=Inf, alpha=0.05, fill=scales::muted("green")) +
   
   scale_colour_manual(values = c("Nucleic acids" = "#8B1914",
                                  "Vitamins and Cofactors" = "#85CDD2"))+ 
   theme_bw()+
   theme(legend.position = "none",
         panel.border = element_rect(colour = "black", fill=NA, size=1),
         text = element_text(size = 25),
         panel.background = element_blank()
   )+ggtitle("Metabolites Perturbed by PRDX1 KO or Etop Treatment")
 ggsave(here::here(output_folder,glue::glue("Metabolite_Behaviour_volcano.pdf")),height = 15,width = 15)
 interesting_metabolites <- c("adenosine_monophosphate","uridine_monophosphate","inosine_triphosphate","guanosine_triphosphate", "aspartic_acid","n_acetylglutamic_acid")
 list_plots_metabo <- list()
 for(i in interesting_metabolites){
   list_plots_metabo[[i]] <-  testing   %>% subset(Metabolite ==  i) %>% 
     mutate(Metabolite = factor(Metabolite,levels = interesting_metabolites),
            Time = case_when(
              Time == "-3" ~ 0,
              Time == "0" ~ 3,
              Time == "24h" ~ 24,
              Time == "8h" ~ 8)) %>% 
     ggplot(aes(x = Time, y = Abundance))   +
     geom_point(alpha = 0.5)+
     #geom_line()+
     theme_bw()+   theme(legend.position = "none",plot.title = element_text(hjust = 0.5),
                         panel.border = element_rect(colour = "black", fill=NA, size=1),
                         text = element_text(size = 35),
                         panel.background = element_blank()
     )+                            # Add confidence intervals
     geom_smooth(aes(fill = Background),linetype=0)+
     scale_fill_manual(values = c(prdx1ko = "#D85E5F",
                                  wt = "#112E88"))+
     ggtitle(i)
 }
 library(gridExtra)
 grid.arrange(grobs = list_plots_metabo, ncol = 2) ## display plot
 ggsave(file = here::here("Output","Figures","Selected_metabolites_aspartic.pdf"),
        arrangeGrob(grobs = list_plots_metabo, ncol = 2),height = 15,width = 15)  ## save plot
 # ggsave(here::here(output_folder,glue::glue("Metabolite_Behaviour_selected.pdf")),height = 15,width = 10)
 list_of_plots <- purrr::imap(.x=DF_anova %>% subset(Significant == T) %>% pull(Metabolite,Significant_type),
                       ~testing %>% subset(Metabolite ==.x) %>% 
                           ggplot(aes(x = Time, y = Abundance, color = Background, 
                                      group = interaction(Background,Replicate)))   +
                           geom_point(alpha = 0.5)+
                           geom_line(alpha = 0.5)+
                           theme_bw()+
                           ggtitle(glue::glue("{.x}")))
 list_of_plots_not_significant <- purrr::imap(.x=DF_anova %>% subset(Significant == F) %>% pull(Metabolite,Significant_type),
                              ~testing %>% subset(Metabolite ==.x) %>% 
                                ggplot(aes(x = Time, y = Abundance, color = Background, 
                                           group = interaction(Background,Replicate)))   +
                                geom_point(alpha = 0.5)+
                                geom_line(alpha = 0.5)+
                                theme_bw()+
                                ggtitle(glue::glue("{.x}")))
 annotate_figure(ggarrange(plotlist  = list_of_plots %>% purrr::keep(.,str_detect(names(.),"and")), ncol = 4,  nrow = 4, common.legend = T),
                 top = text_grob("Significant Treatment and Background Metabolite_Behaviour", color = "black", face = "bold", size = 14))
 ggsave(here::here(output_folder,glue::glue("Significant Treatment and Background Metabolite_Behaviour.pdf")),height = 15,width = 15)
 annotate_figure(ggarrange(plotlist  = list_of_plots %>% purrr::keep(.,names(.) == "Significant Background"), ncol = 4,  nrow = 4, 
                           common.legend = T),                
                 top = text_grob("Significant Background Metabolite_Behaviour", color = "black", face = "bold", size = 14))

 ggsave(here::here(output_folder,glue::glue("Significant Background Metabolite_Behaviour.pdf")),height = 15,width = 15)
 annotate_figure(ggarrange(plotlist  = list_of_plots %>% purrr::keep(.,names(.) == "Significant Treatment"), ncol = 4,  nrow = 3, 
                           common.legend = T), top = text_grob("Significant Treatment Metabolite_Behaviour", color = "black", face = "bold", size = 14))
 ggsave(here::here(output_folder,glue::glue("Significant Treatment Metabolite_Behaviour.pdf")),height = 15,width = 15)
 annotate_figure(ggarrange(plotlist  = list_of_plots_not_significant, ncol = 9,  nrow = 9, common.legend = T),
                 top = text_grob("Non_Significant", color = "black", face = "bold", size = 14))
 ggsave(here::here(output_folder,glue::glue("Significant Treatment and Background Metabolite_Behaviour.pdf")),height = 15,width = 15)
 
 heatmap_df <- Metabolomics%>%# log2() %>% 
    #impute::impute.knn() %>% .[["data"]] %>% 
     as.data.frame() %>% 
    
    rownames_to_column("Metabolite") %>%    rowwise() %>% 
    mutate(Ctrl= base::mean(c(wt_ut_1, wt_ut_2, wt_ut_3), na.rm = T)) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric),~.x-Ctrl))%>% dplyr::select(-Ctrl) %>% 
    pivot_longer(-Metabolite, names_to = "Sample",values_to = "Abundance") %>% 
    mutate(Sample = Sample %>% str_replace_all("prdx1_ko","prdx1ko") %>% 
               str_replace_all("ut_","ut_0_" )) %>% 
    #pivot_longer(contains("_"), names_to = "Sample",values_to = "Abundance") %>% 
    mutate(condition = str_remove_all(Sample,"_.$")) %>% 
    group_by(condition,Metabolite) %>% 
    dplyr::summarise(Abundance_Condition= mean(Abundance)) %>% 
    pivot_wider(names_from = condition, values_from = Abundance_Condition) %>% 
    column_to_rownames("Metabolite") 
    
    #dplyr::select(wt_ut_0,prdx1ko_ut_0,prdx1ko_etop_0,prdx1ko_etop_8h,prdx1ko_etop_24h) #%>%
    
 heatmap_df_etop <- heatmap_df %>%
      rownames_to_column("Metabolite") %>% 
     subset(Metabolite %in% (DF_anova %>% subset(Significant_type == "Significant Treatment") %>% pull(Metabolite)))%>% 
          dplyr::select(Metabolite,wt_etop_0,wt_etop_8h,wt_etop_24h) %>% 
     remove_rownames() %>% 
     column_to_rownames("Metabolite")
     

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(heatmap_df_etop), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(heatmap_df_etop)/paletteLength, max(heatmap_df_etop), length.out=floor(paletteLength/2)))

# Plot the heatmap

pheatmap::pheatmap(heatmap_df_etop,cluster_cols = F,color=myColor, breaks=myBreaks, main = "Metabolites Changing upon etop Treatment and Release" )
##### Diff Corre analysis ####
design_mat_prdx1 = data.frame(prdx1 = c(rep(1,4), rep(0,4)),
                          wt = c(rep(0,4),rep(1,4))) %>% as.matrix()
ddcor_res = ddcorAll(inputMat = heatmap_df %>% dplyr::select(!contains("noco")), design = design_mat_prdx1,
                     compare = c("wt", "prdx1"),corrType = "pearson",dCorAvgMethod = "mean",
                     adjust = "none", nPerm = 0, nPairs = 100)
head(ddcor_res)

heatmap_df%>% dplyr::select(!contains("noco")) %>% rownames_to_column("Metabolite") %>% 
    pivot_longer(-Metabolite, names_to = "Condition", values_to = "Abundance")%>% 
    separate(Condition, c("Background","Treatment","Time")) %>% 
    subset(Metabolite %in% c("deoxycytidine_monophosphate","phosphoenolpyruvic_acid")) %>% 
    mutate(Time = if_else(Treatment == "ut","ut",Time),
           Time = factor(Time, levels = c("ut","0","8h","24h"))) %>% 
    ggplot(aes(x = Time, y = Abundance, color = interaction(Background,Treatment), 
               group = interaction(Background)))   +
    geom_point(alpha = 0.5)+
    geom_line(alpha = 0.5)+
    facet_wrap(~Metabolite)+
    ggtitle("Metabolite Behaviour affected by PRDX1 KO and Etop Treatment")
#### FELLA ####
library(FELLA)
set.seed(1)
#Filter overview pathways
Pathways <- read.delim(here::here("Datasets","Raw","hsa_pathways.txt")) %>% subset(Path_type=="metabolic") %>% 
    pull(Path_id)%>% str_remove_all("hsa")
graph <- buildGraphFromKEGGREST(
        organism = "hsa",
        filter.path = Pathways)

tmpdir <- paste0(tempdir(), "/my_database")
# Mke sure the database does not exist from a former vignette build
# Otherwise the vignette will rise an error
# because FELLA will not overwrite an existing database
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(
     keggdata.graph = graph,
     databaseDir = tmpdir,
     internalDir = FALSE,
     matrices = "pagerank",
     normality = "pagerank",
     niter = 50)
fella.data <- loadKEGGdata(
     databaseDir = tmpdir,
     internalDir = FALSE,
     loadMatrix = "pagerank"
    )
#data("FELLA.sample")
class(fella.data)
#data("input.sample")
#input.full <- c(input.sample, paste0("intruder", 1:10))
compoundsBackground.input <- Mapping$KEGG
input.full = subset(Mapping, Metabolite %in% (DF_anova %>% subset(Significant == T) %>% pull(Metabolite,Significant_type))# list_of_plots$`Significant Background` %>% names() #rownames(annotation_row)[annotation_row$wt_etop_24h_vs_wt_etop_8h == "Significant"]
                    ) %>% pull(KEGG)
show(input.full)
#myAnalysis <- defineCompounds(
#    compounds = input.full, 
#    data = fella.data)
#n.nodes = 300
myAnalysis <- enrich(
    compoundsBackground  = compoundsBackground.input,
    compounds = input.full, 
    method = 'pagerank', 
    approx = "normality", 
    data = fella.data,
    niter=10000)

show(myAnalysis)
plot(
    x = myAnalysis, 
    method = "hypergeom", 
    main = "My first enrichment using the hypergeometric test in FELLA", 
    threshold = 1, 
    data = FELLA.sample)
nlimit  = 150
FELLA::plot(
    x = myAnalysis, 
    method = "pagerank", 
    main = "My first enrichment using the diffusion analysis in FELLA", 
    threshold = 0.05, nlimit = nlimit,
    data = fella.data)
g <- generateResultsGraph(
     object = myAnalysis,
     method = "pagerank",
     nlimit = 300, threshold = 0.05,
     data = fella.data)
plotGraph(
    g)
VisGraph <- visNetwork::toVisNetworkData(g)
single_nodes <- c(VisGraph$edges$from,VisGraph$edges$to) %>% table() %>% 
    .[.<1] %>% names() #%>% str_subset("^R|\\.")
VisGraph$edges <- VisGraph$edges %>% subset(
    !(from %in% single_nodes)|!(to %in% single_nodes)
) %>% mutate(arrows = "to")
VisGraph$nodes <- VisGraph$nodes %>% 
    dplyr::select(-label) %>% 
    dplyr::rename(group = com,
           label = NAME) %>% 
    mutate(
        value= case_when(
            group == 1 ~ 20,
            group == 2 ~ 16,
            group == 3 ~ 4,
            group == 4 ~ 3,
            group == 5 ~ 2,
            TRUE ~ 10),
        color = case_when(
        group == 1 ~ "#BB5151",
            group == 2 ~ "#D2ACD2",
            group == 3 ~ "#E4CB63",
            group == 4 ~ "#B2CEDC",
            group == 5 ~ "#83A483",
            TRUE ~ "ELSE"),
        group = case_when(
        group == 1 ~ "Pathway",
        group == 2 ~ "Module",
        group == 3 ~ "Enzyme",
        group == 4 ~ "Reaction",
        group == 5 ~ "Metabolite",
        TRUE ~ "ELSE"),
        label = label %>% purrr::map_chr(.x = .,~head(.x,1)) %>% str_remove_all(" - Homo sapiens \\(human\\)")) %>%
    subset((!(id %in% single_nodes)))
#    view(VisGraph$nodes)
library(visNetwork)
visNetwork::visNetwork(nodes = VisGraph$nodes, edges = VisGraph$edges, main = "Etoposide WT vs T0") %>% 
    visOptions(highlightNearest = TRUE,selectedBy = "group",collapse = TRUE) %>%
        visLayout(randomSeed = 123)

sheet_names <- openxlsx::getSheetNames(here::here("Datasets","Processed","DEPS_metabolomics.xlsx"))
Comparisons_list <- purrr::map(.x = sheet_names,
                     ~openxlsx::read.xlsx(here::here("Datasets","Processed","DEPS_metabolomics.xlsx"), sheet = .x)) %>% 
  set_names(sheet_names)

Comparisons_list$prdx1_ko_ut_vs_wt_ut_diff %>%
    ggplot(aes(x = log2_FC, y = -log10(p.val), label = Metabolite,  alpha = significant))+
    geom_point(size = 3)+
    ggrepel::geom_text_repel(data = . %>% subset(significant == T),size = 8)+
    #ggrepel::geom_label_repel(data = . %>% subset(significant == T&  Sabatini==T))+
    ggtitle(glue::glue("Diff abundance ", "prdx1_ko_ut_vs_wt_ut"))+
  theme_bw()+theme( text = element_text(size = 25),legend.position = "none")
ggsave(here::here("Output","Figures","Prdx1_ko_vsWT.pdf"),width = 10,height = 10)

###
###KEGG pathway representation####
####Importing Pathways####
pacman::p_load(tidyverse,visNetwork,BiocManager,data.table,
               KEGGgraph, matrixStats,MetaboSignal)
Metabo_Mapping <- fread(here::here("Datasets","Processed","Metabo_Mapping.csv"))
HUMAN_9606_idmapping <- readr::read_tsv(here::here("Datasets","Raw","HUMAN_9606_idmapping.dat"),
                                        col_names = FALSE)
HUMAN_9606_idmapping_hsa <- HUMAN_9606_idmapping[str_detect(HUMAN_9606_idmapping$X3, "hsa:"),-2] %>% 
  left_join(HUMAN_9606_idmapping[HUMAN_9606_idmapping$X2 == "Gene_Name", -2], by = "X1") %>% 
  na.omit() %>% dplyr::select(-X1) %>% set_names(c("KEGG_ID", "Gene_names"))

### Retrieving and Merging KEGG Folate ####
tmp <- tempfile()
Metabo_comparisons <- openxlsx::getSheetNames(here::here("Datasets","Processed","DEPS_metabolomics.xlsx"))
Comparisons_list <- purrr::map(.x = 1:length(Metabo_comparisons),
                               ~openxlsx::read.xlsx(here::here("Datasets","Processed","DEPS_metabolomics.xlsx"),sheet = .x) )%>% 
  set_names(Metabo_comparisons)
pathways <-  c("hsa00240")
#Downloads the genes and compounds in each reaction from KEGG#
Retrieve_genes<- function(pathway){
  retrieveKGML(pathway, organism="hsa", destfile=tmp, method="auto", quiet=TRUE)
  Pathway_genes <- KEGGgraph::parseKGML2Graph(tmp)
  genes_reactions <- data.frame(Gene_id = NULL,
                                Reaction = NULL,
                                pathway = NULL,
                                stringsAsFactors = F)
  for (i in 1:length(Pathway_genes@nodeData@defaults[["KEGGNode"]][["nodes"]])){
    df <- data.frame(Gene_id = Pathway_genes@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@entryID,
                     Reaction = Pathway_genes@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@reaction,
                     pathway = pathway,
                     stringsAsFactors = FALSE)
    genes_reactions <- rbind(genes_reactions,df)
    
  }
  genes_reactions
  
  
}
Retrieve_compounds <- function(pathway){
  retrieveKGML(pathway, organism="hsa", destfile=tmp, method="auto", quiet=TRUE)
  Pathway_metabolites <- KEGGgraph::parseKGML(tmp)
  
  
  compound_reactions <- data.frame(Reaction = NULL,
                                   Substrate = NULL,
                                   Product = NULL,
                                   Direction = NULL,
                                   pathway = NULL,
                                   stringsAsFactors = F)
  
  for (i in 1:length(Pathway_metabolites@reactions)){
    df <- data.frame(Reaction = Pathway_metabolites@reactions[[i]]@name,
                     Substrate = Pathway_metabolites@reactions[[i]]@substrateName,
                     Product = Pathway_metabolites@reactions[[i]]@productName,
                     Direction = Pathway_metabolites@reactions[[i]]@type,
                     pathway = pathway,
                     stringsAsFactors = FALSE)
    compound_reactions <- rbind(compound_reactions,df)
    
  }
  compound_reactions
  
  
}
genes_reactions <- purrr::map_df(pathways,Retrieve_genes)
compound_reactions <-  purrr::map_df(pathways,Retrieve_compounds)

####
Enzymes_Metabolites <- read_tsv(here::here("Datasets","Raw","all_unique_KEGG_metabolites_mapped_KEGG.tsv"))
Metabolites <- Enzymes_Metabolites %>% dplyr::select(KEGG,Compound) %>% distinct
compound_reactions <- separate_rows(compound_reactions, Reaction, sep =" ") %>% 
  mutate(Substrate = str_remove(Substrate,"cpd:"),
         Product = str_remove(Product,"cpd:"))

#Adds Reactions which were missing from KEGGPAth
genes_reactions <- 
  # rbind(genes_reactions, data.frame(Gene_id = c("hsa:1719", "hsa:200895","hsa:4522", "hsa:441024", "hsa:10797","hsa:10797","hsa:4522","hsa:10841","hsa:10841"),
  #                                                  Reaction = c("rn:R00937 rn:R02235", "rn:R00937 rn:R02235","rn:R00943", "rn:R01655", "rn:R01655", "rn:R01218", "rn:R01218","rn:R02302","rn:R03189"),
  #                                                  pathway = "hsa00670")) %>% 
  separate_rows(genes_reactions,Reaction, sep =" ")

Pathways <-     full_join(genes_reactions,compound_reactions, by = c("Reaction", "pathway")) %>% 
  distinct() %>%
  left_join(HUMAN_9606_idmapping_hsa, by= c("Gene_id" = "KEGG_ID")) %>%
  left_join(Metabolites  %>% dplyr::rename(Substrate = KEGG) ) %>% dplyr::rename(From = Compound) %>% 
  left_join(Metabolites  %>% dplyr::rename(Product = KEGG) ) %>% dplyr::rename(To = Compound) %>% 
  distinct(Gene_id, pathway, Substrate, Product, Direction, Gene_names, .keep_all = T) %>%
  mutate(#Gene_name = if_else(is.na(Gene_names),Gene_id,Gene_names),
    From = if_else(is.na(From),Substrate,From),
    To = if_else(is.na(To),Product,To)) %>% 
  mutate(Gene_names = case_when(
    Substrate ==  "C00169" & Product ==  "C00438" ~"CAD",
    Substrate ==  "C00055" & Product ==  "C00112" ~"CMPK1;CMPK2",
    Substrate ==  "C00337" & Product ==  "C00438" ~"CAD",
    Substrate ==  "C01103" & Product ==  "C00105" ~"UMPS",
    Substrate ==  "C00063" & Product ==  "C00055" ~"ENPP1;ENPP3;ASMTL",
    Substrate ==  "C00112" & Product ==  "C00055" & Direction == "irreversible"~"ENTPD8;ENTPD1;ENTPD3",
    Substrate ==  "C00112" & Product ==  "C00055" & Direction == "reversible"~"CMPK1;CMPK2",
    Substrate ==  "C00063" & Product ==  "C00112" & Direction == "irreversible"~"ENTPD8;ENTPD1;ENTPD3",
    Substrate ==  "C00112" & Product ==  "C00063" & Direction == "reversible"~"NME6;AK9;NME7;NME1;NME2;NME3;NME4",
    Substrate ==  "C00475" & Product ==  "C00055" ~"UCKL1;UCK1;UCK2",
    Substrate ==  "C00112" & Product ==  "C00705" ~"RRM2B;RRM1;RRM2",
    Substrate ==  "C00705" & Product ==  "C00458" ~"NME6;AK9;NME7;NME1;NME2;NME3;NME4",
    Substrate ==  "C00239" & Product ==  "C00705" ~"CMPK2;CMPK1",
    Substrate ==  "C00239" & Product ==  "C00881" ~"HDDC2;NT5C2;NT5DC4;NT5C;NT5E;NT5M",
    Substrate ==  "C00881" & Product ==  "C00526" ~"CDA",
    Substrate ==  "C00365" & Product ==  "C00526" ~"TK1;TK2",
    Substrate ==  "C00429" & Product ==  "C00106" ~"DPYD",
    Substrate ==  "C00429" & Product ==  "C02642" ~"DPYS",
    Substrate ==  "C02642" & Product ==  "C00099" ~"UPB1",
    Substrate ==  "C01346" & Product ==  "C00460" ~"NME6;AK9;NME7;NME1;NME2;NME3;NME4",
    Substrate ==  "C00363" & Product ==  "C00364" ~"ENTPD8;ENTPD1;ENTPD3",
    Substrate ==  "C00364" & Product ==  "C00363" ~"DTYMK",
    Substrate ==  "C00214" & Product ==  "C00364" ~"TK1;TK2",
    Substrate ==  "C00459" & Product ==  "C00364" ~"ENTPP1;ENTPP3l;ASMTL",
    Substrate ==  "C00363" & Product ==  "C00459" ~"NME6;AK9;NME7;NME1;NME2;NME3;NME4",
    Substrate ==  "C00459" & Product ==  "C00363" ~"ENTPD8;ENTPD1;ENTPD3",
    
    Substrate ==  "C00214" & Product ==  "C00672" ~"TYMP",
    Substrate ==  "C00214" & Product ==  "C00178" ~"TYMP",
    TRUE ~ Gene_names
  )) %>% separate_rows(Gene_names,sep = ";")


####Network database Creating ####
nodes <- data.frame(id =c(Pathways$From,Pathways$To),
                    label =c(Pathways$From,Pathways$To),
                    KEGG= c(Pathways$Substrate,Pathways$Product)) %>% 
  subset(!duplicated(id))
# mypal <- colorRampPalette( c( "#053061"   ,"#F7F7F7","#670A1F") )#(10)
vector_low_mid_high <- c( "#053061"   ,"#F7F7F7","#670A1F")
# map2color<-function(x,pal,limits=NULL){
#     if(is.null(limits)) limits=range(x)
#     pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
# }
map2color_center<-function(x,vector_low_mid_high){
  # if(is.null(limits)) limits=range(x)
  # vector_low_mid_high <- c(low = "blue",
  #                          mid = "white",
  #                          high = "red")
  # pal = mypal
  # x = Volcano_DFs$x0h_vs_x24h_diff$log2_FC
  # vector_low_mid_high = c("blue","white","red")
  # x = Volcano_DFs$dmso_vs_x0h_diff$log2_FC 
  names(x) <- 1:length(x)
  pal_high <- colorRampPalette( vector_low_mid_high[-1] )(8 )
  pal_low <- colorRampPalette( vector_low_mid_high[-3] )(8 )
  x_high <- x[x>0]
  x_low <- x[x<=0]
  limits_high=range(x_high)
  limits_low=range(x_low)
  x_high <-  pal_high[findInterval(x_high,seq(limits_high[1],limits_high[2],length.out=length(pal_high)+1), all.inside=TRUE)] %>% 
    set_names(names(x_high))
  x_low <-  pal_low[findInterval(x_low,seq(limits_low[1],limits_low[2],length.out=length(pal_low)+1), all.inside=TRUE)]%>% 
    set_names(names(x_low))
  c(x_high,x_low) %>% .[names(x)] %>% as.vector() 
  
}

# Volcano_DFs <- map(.x = 1:3,
#                    ~openxlsx::read.xlsx(here::here("Datasets","Processed", "Batch_corrected_DEPs.xlsx"), sheet = .x, rowNames = T) %>% 
#                        subset(!is.na(ID))) %>% 
#     set_names(openxlsx::getSheetNames(here::here("Datasets","Processed", "Batch_corrected_DEPs.xlsx")))

edges_enzymes_omics <- purrr::map(.x =Volcano_DFs,
                                  ~
                                    # .x = Volcano_DFs[[1]]
                                    Pathways %>% 
                                    dplyr::select(c( "From","To", "Gene_names", "Direction"))%>% 
                                    dplyr::rename(label = Gene_names, from = From, to = To) %>% 
                                    left_join(.x %>% mutate(
                                      color = map2color(log2_FC,mypal)),
                                      
                                      , by = c("label" = "ID")) %>% 
                                    mutate(font.color = "green",
                                           
                                           font.size = 10,
                                           width = if_else(is.na(color),2,5),
                                           # label = if_else(is.na(label)," ",label),
                                           color = if_else(is.na(log2_FC),"#a6a6a6", color),
                                           
                                           arrows.middle = if_else(Direction == "irreversible", T, F),
                                           dashes = if_else(label == " ", T, F) ) %>%
                                    distinct(from,to,label,.keep_all = T)) %>% set_names(names(Volcano_DFs))
nodes_omics <- purrr::map(.x =Comparisons_list,
                          ~
                            
                            nodes %>% 
                            left_join(.x %>% mutate(
                              color = map2color(log2_FC,mypal)),
                              
                              , by = c("KEGG")) %>% 
                            mutate(
                              font.color = "black",
                              font.size = 10,
                              color = if_else(is.na(log2_FC),"#E7E7E7", color)
                              
                              # arrows.middle = if_else(Direction == "irreversible", T, F),
                              # dashes = if_else(label == " ", T, F)
                            ))  %>% set_names(names(Comparisons_list))

#nodes <- nodes %>% mutate(id = MetaboSignal::MS_changeNames(paste0("cpd:", id),"hsa"))
list_nodes <- list(1,1,2,2,3,3,4,4)
list_edges <- list(rep(c(4,5),times = 4)) %>% unlist()
purrr::walk2(list_nodes,list_edges, make_visnetwork)
make_visnetwork <- function(x,y){
  print(x)
  print(y)
  # behaviour <- if_else(y == 4,"")
  graph_test <- visNetwork(nodes_omics[[x]],edges_enzymes_omics[[y]], 
                           height = 1000, width = "100%",
                           main = paste("Metabolic Vulnerabolities - ",names(edges_enzymes_omics)[y],names(nodes_omics)[x]),
                           submain ="Essential enzymes for clearance have high sgRNA beta (red), Essential metabolites are upregulated (red) <br>Red metabolite = Upreg <br>Red enzyme = more guides in treatment = enzyme necessary for clearance?") %>% 
    visIgraphLayout(smooth = T,randomSeed = 123)%>% 
    visNodes(
      shape = "box",
      shadow = list(enabled = TRUE, size = 10)
    ) %>%
    visEdges(arrows = list(middle = list(scaleFactor = 0.1))) %>% 
    visEdges(smooth = list(type = "dynamic"), color = list(highlight = "#C62F4B",opacity = 0.35, border = "black"), 
             font = list(align = "middle")) %>%
    visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
               selectedBy = "KEGG") %>%
    visPhysics(stabilization = FALSE,solver = "forceAtlas2Based", 
               forceAtlas2Based = list(gravitationalConstant = -150)) 
  visSave(graph_test, here::here("Output",glue::glue("{names(edges_enzymes_omics)[y]},{names(nodes_omics)[x]},.html")), 
          selfcontained = TRUE, background = "white")}
pathways <- c("hsa00020","map00250")
Retrieve_genes<- function(pathway){
  retrieveKGML(pathway, organism="hsa", destfile=tmp, method="auto", quiet=TRUE)
  Pathway_genes <- KEGGgraph::parseKGML2Graph(tmp)
  genes_reactions <- data.frame(Gene_id = NULL,
                                Reaction = NULL,
                                pathway = NULL,
                                stringsAsFactors = F)
  for (i in 1:length(Pathway_genes@nodeData@defaults[["KEGGNode"]][["nodes"]])){
    df <- data.frame(Gene_id = Pathway_genes@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@entryID,
                     Reaction = Pathway_genes@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@reaction,
                     pathway = pathway,
                     stringsAsFactors = FALSE)
    genes_reactions <- rbind(genes_reactions,df)
    
  }
  genes_reactions
  
  
}
Retrieve_compounds <- function(pathway){
  retrieveKGML(pathway, organism="hsa", destfile=tmp, method="auto", quiet=TRUE)
  Pathway_metabolites <- KEGGgraph::parseKGML(tmp)
  
  
  compound_reactions <- data.frame(Reaction = NULL,
                                   Substrate = NULL,
                                   Product = NULL,
                                   Direction = NULL,
                                   pathway = NULL,
                                   stringsAsFactors = F)
  
  for (i in 1:length(Pathway_metabolites@reactions)){
    df <- data.frame(Reaction = Pathway_metabolites@reactions[[i]]@name,
                     Substrate = Pathway_metabolites@reactions[[i]]@substrateName,
                     Product = Pathway_metabolites@reactions[[i]]@productName,
                     Direction = Pathway_metabolites@reactions[[i]]@type,
                     pathway = pathway,
                     stringsAsFactors = FALSE)
    compound_reactions <- rbind(compound_reactions,df)
    
  }
  compound_reactions
  
  
}
genes_reactions <- purrr::map_df(pathways,Retrieve_genes)
compound_reactions <-  purrr::map_df(pathways,Retrieve_compounds)

####
Enzymes_Metabolites <- read_tsv(here::here("Datasets","Raw","all_unique_KEGG_metabolites_mapped_KEGG.tsv"))
Metabolites <- Enzymes_Metabolites %>% dplyr::select(KEGG,Compound) %>% distinct
compound_reactions <- separate_rows(compound_reactions, Reaction, sep =" ") %>% 
  mutate(Substrate = str_remove(Substrate,"cpd:"),
         Product = str_remove(Product,"cpd:"))

#Adds Reactions which were missing from KEGGPAth
genes_reactions <- 
  # rbind(genes_reactions, data.frame(Gene_id = c("hsa:1719", "hsa:200895","hsa:4522", "hsa:441024", "hsa:10797","hsa:10797","hsa:4522","hsa:10841","hsa:10841"),
  #                                                  Reaction = c("rn:R00937 rn:R02235", "rn:R00937 rn:R02235","rn:R00943", "rn:R01655", "rn:R01655", "rn:R01218", "rn:R01218","rn:R02302","rn:R03189"),
  #                                                  pathway = "hsa00670")) %>% 
  separate_rows(genes_reactions,Reaction, sep =" ")

Pathways <-     full_join(genes_reactions,compound_reactions, by = c("Reaction", "pathway")) %>% 
  distinct() %>%
  left_join(HUMAN_9606_idmapping_hsa, by= c("Gene_id" = "KEGG_ID")) %>%
  left_join(Metabolites  %>% dplyr::rename(Substrate = KEGG) ) %>% dplyr::rename(From = Compound) %>% 
  left_join(Metabolites  %>% dplyr::rename(Product = KEGG) ) %>% dplyr::rename(To = Compound) %>% 
  distinct(Gene_id, pathway, Substrate, Product, Direction, Gene_names, .keep_all = T) %>%
  mutate(#Gene_name = if_else(is.na(Gene_names),Gene_id,Gene_names),
    From = if_else(is.na(From),Substrate,From),
    To = if_else(is.na(To),Product,To))

New_map <- list(edges= data.frame(from = c("C00122","C00122","C00149","C00036",
                                           "C03406","C03406","C00062","C00077","C00327",
                                           "C00042","C00122_m","C00149_m","C00036_m","C00036_m",
                                           "C00158","C00311","C00026","C00042",
                                           "C00025","C00064","C00025","C00624",
                                           "C00169","C00169","C00169","C00049","C00049","C00049","C00036_m","C00149_m",
                                           "C00026","C00064"),
                                  to = c("C00122_m","C00149","C00149_m","C00149",
                                         "C00122","C00062","C00077","C00327","C03406",
                                         "C00122_m","C00149_m","C00036_m","C03406","C00158",
                                         "C00311","C00026","C00042","C00122_m",
                                         "C00026","C00025","C00624","C00169",
                                         "C00438","C00327","C00438","C00438","C03406","C00036","C00049","C00036_m",
                                         "C00042","C00064_m")) %>% 
                  distinct(),
                nodes = data.frame(id = c("C00122","C00122_m","C00149","C00149_m","C00036","C00036_m","C03406",
                                          "C00062","C00049","C00327","C00077","C00158",
                                          "C00311","C00042","C00026","C00025","C00169",
                                          "C00624","C00438","nucleotides","C00064_m","C00064"),
                                   label = c("Fumarate","Fumarate_m","Malate","Malate_m","Oxaloacetate","Oxaloacetate_m","Argininosuccinic acid",
                                             "Arginine","Aspartate","Citrulline","Ornithine","Citrate",
                                             "Isocitrate","Succinate","alpha-Ketoglutaric acid","Glutamate","Carbamoyl phosphate",
                                             "N-Acetyl-L-glutamate","N-Carbamoyl-L-aspartate","nucleotides","Glutamine_m","Glutamine")))
New_map$nodes <- New_map$nodes %>% mutate(KEGG= str_remove_all(id,"_m"))
New_map$edges <- New_map$edges %>% mutate(Substrate =str_remove_all(from,"_m") ,
                                          Product =str_remove_all(to,"_m") )
edges <- left_join(New_map$edges,Pathways %>% dplyr::select(Substrate,Product,Gene_names,Direction)   )%>% 
  mutate(Gene_names = case_when(
    Substrate ==  "C00122" & Product ==  "C00149" ~"FH",
    from ==  "C00149" & to ==  "C00149_m"~"SLC25A11;SLC25A1",
    # from ==  "C00149" & to ==  "C00149_m"~"SLC25A1",
    from ==  "C00064" & to ==  "C00064_m"~"SLC1A5",
    from ==  "C00064" & to ==  "C00025"~"GLS;GLS2",
    from ==  "C00025" & to ==  "C00026"~"GLUD1;GLUD2",
    Substrate ==  "C00158" & Product ==  "C00311" ~"ACO1;ACO2",
    Substrate ==  "C00169" & Product ==  "C00327" ~"OTC",
    Substrate ==  "C00077" & Product ==  "C00327" ~"OTC",
    Substrate ==  "C00169" & Product ==  "C00438" ~"CAD",
    Substrate ==  "C03406" & Product ==  "C00062" ~"ASL",
    Substrate ==  "C00327" & Product ==  "C03406" ~"ASS1",
    Substrate ==  "C00062" & Product ==  "C00077"~"ARG1;ARG2",
    Substrate ==  "C00025" & Product ==  "C00624" ~"NAGS",
    Substrate ==  "C00624" & Product ==  "C00169" ~"CPS1",
    Substrate ==  "C00149" & Product ==  "C00036" ~"MDH1;MDH2",
    Substrate ==  "C00026" & Product == "C00042"~"OGDH;DLST;SUCLG2;SUCLG1",
    TRUE ~ Gene_names
  )) %>% separate_rows(Gene_names,sep = ";")
screen_list <- list(D10 =  Screen_input_volc %>% dplyr::select(Gene,`high|beta`) %>% dplyr::rename(ID = Gene,log2_FC = `high|beta`),
                    D14 = norm_genes%>% dplyr::select(Gene,Etop_vs_DMS0)%>% dplyr::rename(ID = Gene,log2_FC = Etop_vs_DMS0))
edges <- purrr::map(.x = screen_list,
                    ~
                      # .x = Volcano_DFs[[1]]
                      edges %>% 
                      # dplyr::select(c( "From","To", "Gene_names", "Direction"))%>% 
                      dplyr::rename(label = Gene_names) %>% 
                      left_join(.x %>% subset(!is.na(ID)) %>%  
                                  mutate(
                                    color = map2color_center(-log2_FC,vector_low_mid_high),
                                    log2_FC_rev = -log2_FC),
                                
                                , by = c("label" = "ID")) %>% 
                      mutate(font.color = "green",
                             
                             font.size = 10,
                             width = if_else(is.na(color),2,5),
                             # label = if_else(is.na(label)," ",label),
                             color = if_else(is.na(log2_FC),"#a6a6a6", color),
                             
                             # arrows.middle = if_else(Direction == "irreversible", T, F),
                             dashes = if_else(label == " ", T, F) ) %>%
                      distinct(from,to,label,.keep_all = T)) %>% set_names(names(screen_list))


nodes <- purrr::map(.x =Comparisons_list, ~
                      
                      New_map$nodes %>% 
                      left_join(.x %>% mutate(
                        color = map2color_center(log2_FC,vector_low_mid_high )),
                        
                        , by = c("KEGG")) %>% 
                      mutate(
                        font.color = "black",
                        font.size = 10,
                        color = if_else(is.na(log2_FC),"#E7E7E7", color)
                        
                        # arrows.middle = if_else(Direction == "irreversible", T, F),
                        # dashes = if_else(label == " ", T, F)
                      ))  %>% set_names(names(Comparisons_list))
visNetwork(nodes =  nodes$wt_etop_24h_vs_wt_ut_diff ,edges = edges$D14, 
           height = 1000, width = "100%", 
           main = "Metabolic changes 24hr vs DMSO, essentiality D5 recovery vs DMSO") %>% 
  visEdges(arrows = list(middle = list(scaleFactor = 1)))


