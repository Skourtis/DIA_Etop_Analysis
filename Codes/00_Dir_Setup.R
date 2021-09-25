#Settting up directory
##Settting up directory
options(repos = getOption("repos")["CRAN"])
# install.packages("PTXQC")
install.packages("pacman")
pacman::p_load(piggyback, renv, here, tidyverse, targets,
               visNetwork,matrixStats,magick,testthat, openxlsx, janitor,seqinr)
renv::init()
#usethis::use_test()
###Creates all files needed for project
folders <- c("Datasets","Output","Datasets/Raw","Datasets/Processed")
purrr::walk(.x = folders,~dir.create(here::here(.x)))
#piggyback::pb_new_release()

## Created a first release directly on Github
#pb_new_release("Skourtis/Project_Template")
piggyback::pb_track(c("Datasets/Raw/*.txt",
                      "Datasets/Raw/*.dat",
                      "Datasets/Raw/*.zip",
                      "Datasets/Raw/*.csv",
                      "Datasets/Raw/*.RData"))%>%
    pb_upload(repo = "Skourtis/DIA_Etop_Analysis")

##end
renv::snapshot()


