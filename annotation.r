library(data.table)
library(dplyr)


setwd("~/Dropbox/SNU_medical/EWAS") # Linux
# setwd("C:/Users/kevin/Dropbox/SNU_medical/EWAS")
bmiq_data_2yr <- fread("cg_data/2yr/bmiq.norm.data.txt", header = T)
bmiq_data_6yr <- fread("cg_data/6yr/bmiq.norm.data.txt", header = T)

annote_450k <- fread("annotation/annotation_450k.csv", header = T)
annote_epic <- fread("annotation/annotation_epic.csv", header = T)

str(annote_450k)
ess_annote <- annote_epic %>% select(Name, UCSC_RefGene_Name, UCSC_RefGene_Group, Relation_to_UCSC_CpG_Island, SNP_DISTANCE, SNP_MinorAlleleFrequency, CHR)
ess_annote <- annote_450k %>% select(Name, UCSC_RefGene_Name, UCSC_RefGene_Group, Relation_to_UCSC_CpG_Island, SNP_DISTANCE, SNP_MinorAlleleFrequency, CHR)
spsp <- function(x){
  sapply(x %>% strsplit(";"), function(y) y[1])
}

ess_annote1 <- ess_annote %>% transmute(cg_site = Name,
                                        RefGene_Name = spsp(UCSC_RefGene_Name),
                                        RefGene_Group = spsp(UCSC_RefGene_Group),
                                        RefGene_Island = spsp(Relation_to_UCSC_CpG_Island),
                                        SNP_DISTANCE = spsp(SNP_DISTANCE),
                                        SNP_MinorAlleleFrequency = spsp(SNP_MinorAlleleFrequency),
                                        chrome = CHR
)

data_m <- merge(data, ess_annote1, by = "cg_site", all.x = T)

data_m <- data_m %>% select(-V1, -rownumber)

data_arr <- arrange(data_m, p_value)

data_filt <- data_arr %>% subset(!(SNP_DISTANCE %in% c(0,1))) %>%
  subset(SNP_MinorAlleleFrequency <= 0.05) %>%
  subset(!(chrome %in% c("X", "Y")))

data_f <- data_filt %>% mutate(rownumbr = 1:nrow(data_filt),
                               fdr5 = rownumbr*0.05/nrow(data_filt),
                               fdr_result = ifelse(p_value < fdr5, "Y", "N"))



annote_filter <- function(data, annote_epic){
  ess_annote <- annote_epic %>% select(Name, UCSC_RefGene_Name, UCSC_RefGene_Group, Relation_to_UCSC_CpG_Island, SNP_DISTANCE, SNP_MinorAlleleFrequency, CHR)
  
  spsp <- function(x){
    sapply(x %>% strsplit(";"), function(y) y[1])
  }
  
  ess_annote1 <- ess_annote %>% transmute(cg_site = Name,
                                          RefGene_Name = spsp(UCSC_RefGene_Name),
                                          RefGene_Group = spsp(UCSC_RefGene_Group),
                                          RefGene_Island = spsp(Relation_to_UCSC_CpG_Island),
                                          SNP_DISTANCE = spsp(SNP_DISTANCE),
                                          SNP_MinorAlleleFrequency = spsp(SNP_MinorAlleleFrequency),
                                          chrome = CHR
  )
  
  data_m <- data %>% left_join(ess_annote1, by= "cg_site")
  # data_m <- merge(data, ess_annote1, by = "cg_site", all.x = T)
  
  data_m <- data_m %>% select(-rownumber)
  
  data_arr <- arrange(data_m, p_value)
  
  data_filt <- data_arr %>% subset(!(SNP_DISTANCE %in% c(0,1))) %>%
    subset(SNP_MinorAlleleFrequency <= 0.05) %>%
    subset(!(chrome %in% c("X", "Y")))
  
  data_f <- data_filt %>% mutate(rownumbr = 1:nrow(data_filt),
                                 fdr5 = rownumbr*0.05/nrow(data_filt),
                                 fdr_result = ifelse(p_value < fdr5, "Y", "N"))
  return(data_f)
}
