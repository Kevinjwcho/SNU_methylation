# install.packages("BiocManager")
library(BiocManager)
# BiocManager::install("sva")
# BiocManager::install("bladderbatch")
# BiocManager::install("pamr")
# BiocManager::install("limma")
# install.packages("data.table")
# install.packages("dplyr")

library(sva)
# library(bladderbatch)
# library(pamr)
# library(limma)

library(data.table)
library(dplyr)


setwd("C:/Users/kevin/Dropbox/SNU_medical")

## Loading data ----------------
raw_2yr <- fread("EWAS/cg_data/2yr/bmiq.norm.data.txt", header = T)
raw_6yr <- fread("EWAS/cg_data/6yr/bmiq.norm.data.txt", header = T)

raw_2yr_new <- raw_2yr %>% select(-V1) %>% as.matrix
raw_6yr_new <- raw_6yr %>% select(-V1) %>% as.matrix

rownames(raw_2yr_new) <- raw_2yr$V1
rownames(raw_6yr_new) <- raw_6yr$V1

## Sample data ---------------
sam_2yr <- fread("EWAS/cg_data/2yr/Sample.Table.txt", header = T)
sam_6yr <- fread("EWAS/cg_data/6yr/Sample.Table.txt", header = T)

## Annotation --------------

# annote_450k <- fread("annotation/annotation_450k.csv", header = T)
annote_epic <- fread("EWAS/annotation/annotation_epic.csv", header = T)

str(annote_epic)
ess_annote <- annote_epic %>% select(Name, UCSC_RefGene_Name, UCSC_RefGene_Group, Relation_to_UCSC_CpG_Island, SNP_DISTANCE, SNP_MinorAlleleFrequency, CHR)
# ess_annote <- annote_450k %>% select(Name, UCSC_RefGene_Name, UCSC_RefGene_Group, Relation_to_UCSC_CpG_Island, SNP_DISTANCE, SNP_MinorAlleleFrequency, CHR)
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

raw_2yr <- raw_2yr %>% mutate(cg_site = V1) %>% select(-V1)
raw_6yr <- raw_6yr %>% mutate(cg_site = V1) %>% select(-V1)

merge_dat_2yr <- raw_2yr %>% left_join(ess_annote1, by = "cg_site")
merge_dat_6yr <- raw_6yr %>% left_join(ess_annote1, by = "cg_site")

scaling <- raw_2yr %>% select(cg_site) %>% mutate(age2yr = 1) %>%
  inner_join(raw_6yr %>% select(cg_site) %>% mutate(age6yr = 1))

str(scaling)
dim(scaling)
dim(scaling %>% na.omit)

merge_2yr_250k <- merge_dat_2yr %>% filter(cg_site %in% scaling$cg_site) 
merge_6yr_250k <- merge_dat_6yr %>% filter(cg_site %in% scaling$cg_site) 

del_candi <- names(merge_6yr_250k)[82:87]
data_filt_2yr_500k <- merge_dat_2yr %>% subset(!(SNP_DISTANCE %in% c(0,1))) %>%
  subset(SNP_MinorAlleleFrequency <= 0.05) %>% 
  subset(!(chrome %in% c("X", "Y"))) %>% select(-del_candi)

data_filt_2yr_250k <- merge_2yr_250k %>% subset(!(SNP_DISTANCE %in% c(0,1))) %>%
  subset(SNP_MinorAlleleFrequency <= 0.05) %>%
  subset(!(chrome %in% c("X", "Y"))) %>% select(-del_candi)

data_filt_6yr_250k <- merge_6yr_250k %>% subset(!(SNP_DISTANCE %in% c(0,1))) %>%
  subset(SNP_MinorAlleleFrequency <= 0.05) %>%
  subset(!(chrome %in% c("X", "Y"))) %>% select(-del_candi)

## batch effect -------------

batch_2yr <- sam_2yr$`Sentrix Barcode` %>% as.numeric
batch_6yr <- sam_6yr$`Sentrix Barcode` %>% as.numeric
  
# 1로만 이루어진 디자인 매트릭스를 만듭니다.
modcombat_2yr <- model.matrix(~1, data=batch_2yr %>% as.data.frame)
modcombat_6yr <- model.matrix(~1, data=batch_6yr %>% as.data.frame)

# matrix로 데이터를 바꾸기 ---------
matrix_2yr_500k <- data_filt_2yr_500k %>% select(-cg_site) %>% as.matrix
rownames(matrix_2yr_500k) <- data_filt_2yr_500k$cg_site

matrix_2yr_250k <- data_filt_2yr_250k %>% select(-cg_site) %>% as.matrix
rownames(matrix_2yr_250k) <- data_filt_2yr_250k$cg_site

matrix_6yr_250k <- data_filt_6yr_250k %>% select(-cg_site) %>% as.matrix
rownames(matrix_6yr_250k) <- data_filt_6yr_250k$cg_site


# Batch effect correction
combat_2yr_500k <- ComBat(dat = matrix_2yr_500k, batch = batch_2yr, mod = modcombat_2yr, par.prior = T)
combat_2yr_250k <- ComBat(dat = matrix_2yr_250k, batch = batch_2yr, mod = modcombat_2yr, par.prior = T)
combat_6yr_250k <- ComBat(dat = matrix_6yr_250k, batch = batch_6yr, mod = modcombat_6yr, par.prior = T)


# Convert the data into data.frame form.
combat_2yr_500k <- combat_2yr_500k %>% as.data.frame
combat_2yr_250k <- combat_2yr_250k %>% as.data.frame
combat_6yr_250k <- combat_6yr_250k %>% as.data.frame

# 이걸 저장합니다.
fwrite(combat_2yr_500k, "EWAS/Batcheffect/kslee_result/batch_2yr_500k.csv", row.names = TRUE)
fwrite(combat_2yr_250k, "EWAS/Batcheffect/kslee_result/batch_2yr_250k.csv", row.names = TRUE)
fwrite(combat_6yr_250k, "EWAS/Batcheffect/kslee_result/batch_6yr_250k.csv", row.names = TRUE)



