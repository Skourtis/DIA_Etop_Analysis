#Settting up directory
##Settting up directory
# options(repos = getOption("repos")["CRAN"])
# install.packages("PTXQC")
#install.packages("pacman")
# if (!requireNamespace("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
# remotes::install_github("vdemichev/diann-rpackage")
# installr::updateR()
# BiocManager::install("DEP")
# remotes::install_github("demar01/PeCorA")
if(!"RCy3" %in% installed.packages()){
    install.packages("BiocManager")
    BiocManager::install("RCy3")
}
library(RCy3)
library(igraph)
pacman::p_load(piggyback, renv, here, tidyverse, targets, DEP,pheatmap,diann,PeCorA,sva,imp4p,
               org.Hs.eg.db,clusterProfiler,ggridges,usethis,
               visNetwork,matrixStats,magick,testthat, openxlsx, janitor,seqinr)
renv::init()
#usethis::use_test()
###Creates all files needed for project
folders <- c("Datasets","Output","Datasets/Raw","Datasets/Processed")
purrr::walk(.x = folders,~dir.create(here::here(.x)))
#piggyback::pb_new_release()

## Created a first release directly on Github
#pb_new_release("Skourtis/Project_Template")

usethis::use_git_ignore(c("Datasets/Raw/*.txt",
                      "Datasets/Raw/*.dat",
                      "Datasets/Raw/*.tsv",
                      "Datasets/Raw/*.zip",
                      "Datasets/Raw/*.csv",
                      "Output/Etop_DIA_EC_two_methods/*.png",
                      "Datasets/Raw/*.RData"))
    pb_upload(repo = "Skourtis/DIA_Etop_Analysis")

##end
renv::snapshot()


