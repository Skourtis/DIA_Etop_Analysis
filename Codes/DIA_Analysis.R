#### Etop DIA ####
#### Loading####

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
output_folder = here::here("Output","Etop_DIA_EC_106M955")
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
                      Condition = c("DMSO_1","DMSO_2","DMSO_3","0H_1", "0H_2","0H_3","24H_1", "24H_2" ,"24H_3")) %>% 
    mutate(Condition = janitor::make_clean_names(Condition))

Method1_P11833_41 <- Load_DIA_NN_Data(DIA_report_file,Samples)
