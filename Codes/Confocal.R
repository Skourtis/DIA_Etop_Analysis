##### Confocal

COX4 <- openxlsx::read.xlsx(here::here("Datasets","Raw","COX4 - 20x confocal (1).xlsx"))
COX4<- COX4 %>% pivot_longer(everything(),names_to = "Sample",values_to = "Intensity") %>% 
    separate(col = "Sample",into = c("Gene","Time","Subce","Treatment"))
COX4 %>% ggplot(aes(x = interaction(Time,Treatment), y = Intensity))+
    geom_boxplot()+
    facet_wrap("Subce")
