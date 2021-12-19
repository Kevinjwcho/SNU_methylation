library(data.table)
library(dplyr)

setwd("~/Dropbox/SNU_medical")
# setwd("C:/Users/kevin/Dropbox/SNU_medical")

bmiq_data_2yr <- fread("EWAS/MDS_EWAS/batch_2yr_403k.csv", header = T)
bmiq_data_6yr <- fread("EWAS/MDS_EWAS/batch_6yr_403k.csv", header = T)

sample_id <- colnames(bmiq_data_2yr)[-1]

result_dat <- fread("EWAS/MDS_EWAS/report_0605/result.csv", header = T)

result_annote <- result_dat %>%
  transmute(cg_site = cg_site,
            year = year,
            cg_chrome = CHR,
            cg_position = MAPINFO,
            ) %>% na.omit


gene_dat <- fread("gene/genome/20190114_SNUCM_HongYunChul_Genotype.txt/20190114_SNUCM_HongYunChul_Genotype.txt")

gene_id <- colnames(gene_dat)[-c(1:9)]
inter_id <- intersect(sample_id, gene_id)
gene_dat_annote <- gene_dat %>% mutate(SNP = `dbSNP RS ID`,
                                       chrome = Chromosome,
                                       position = `Physical Position` %>% as.numeric) %>% 
  select(SNP, chrome, position, inter_id)


colnames(bmiq_data_2yr)[1] <- "cg_site"
colnames(bmiq_data_6yr)[1] <- "cg_site"

bmiq_data_2yr_sel <- bmiq_data_2yr %>% select(cg_site, inter_id)
bmiq_data_6yr_sel <- bmiq_data_6yr %>% select(cg_site, inter_id)

result <- data.frame()
n <- nrow(result_annote)
for(i in 1:n){
  cg_n <- result_annote[i, 1] %>% as.character()
  if(result_annote[i, 4] == 2){
    sub_methyl <- bmiq_data_2yr_sel %>% filter(cg_site == cg_n) %>% select(-cg_site) %>% t() %>% as.data.frame
    sub_methyl <- sub_methyl %>% mutate(regist_no = rownames(sub_methyl))
  }else{
    sub_methyl <- bmiq_data_6yr_sel %>% filter(cg_site == cg_n) %>% select(-cg_site) %>% t() %>% as.data.frame
    sub_methyl <- sub_methyl %>% mutate(regist_no = rownames(sub_methyl))
  }
  colnames(sub_methyl)[1] <- "methyl"
  sub_result <- gene_dat_annote %>% filter(chrome == result_annote[i, 3] %>% as.numeric,
                                           position <= result_annote[i, 4] %>% as.numeric + 10^5,
                                           position >= result_annote[i, 4] %>% as.numeric - 10^5) %>% 
    transmute(SNP = SNP, SNP_chrome = chrome, SNP_position = position)
  
  sub_result_new <- data.frame(result_annote[i,], sub_result)
  
  
  allele_find <- function(x){
    chrome1<- x[6] %>% as.numeric
    position1<- x[7] %>% as.numeric
    allele_dat <- gene_dat_annote %>% filter(chrome == chrome1, position == position1)
    allele_dat1 <- allele_dat[,-c(1:3)] %>% unlist()
    allele_dat2 <- allele_dat1[which(allele_dat1 != "NN")]
    allele_tab <- table(allele_dat2)
    end <- paste(names(allele_tab), collapse = ", ")
    return(end)
  }
  
  allele_num <- function(x){
    chrome1<- x[6] %>% as.numeric
    position1<- x[7] %>% as.numeric
    allele_dat <- gene_dat_annote %>% filter(chrome == chrome1, position == position1)
    allele_dat1 <- allele_dat[,-c(1:3)] %>% unlist()
    allele_dat2 <- allele_dat1[which(allele_dat1 != "NN")]
    allele_tab <- table(allele_dat2)
    end <- length(allele_tab)
    return(end)
  }
  
  n_obs <- function(x){
    chrome1<- x[6] %>% as.numeric
    position1<- x[7] %>% as.numeric
    allele_dat <- gene_dat_annote %>% filter(chrome == chrome1, position == position1)
    allele_dat1 <- allele_dat[,-c(1:3)] %>% unlist()
    allele_dat2 <- allele_dat1[which(allele_dat1 != "NN")]
    # allele_tab <- table(allele_dat2)
    end <- length(allele_dat2)
    return(end)
  }
  
  sub_result_new$allele <- apply(sub_result_new, 1, allele_find)
  sub_result_new$allele_num <- apply(sub_result_new, 1, allele_num)
  sub_result_new$n_obs <- apply(sub_result_new, 1, n_obs)
  
  
  anova_p <- function(x){
    chrome1<- x[6] %>% as.numeric
    position1<- x[7] %>% as.numeric
    allele_num1<- x[9] %>% as.numeric
    if(allele_num1 == 1){
      return(NA)
    }else{
      allele_dat <- gene_dat_annote %>% filter(chrome == chrome1, position == position1)
      allele_dat1 <- allele_dat[,-c(1:3)] %>% t() 
      allele_dat2 <- allele_dat1 %>% as.data.frame %>% mutate(regist_no = rownames(allele_dat1))
      colnames(allele_dat2)[1] <- "allele"
      sub_dat <- sub_methyl %>% left_join(allele_dat2, by = "regist_no") %>% filter(allele != "NN")
      anov_smry <- lm(methyl ~ allele, data = sub_dat) %>% anova()
      p_v<- anov_smry$`Pr(>F)`[1]
      return(p_v)
    }
  }
  
  kruskal_p <- function(x){
    chrome1<- x[6] %>% as.numeric
    position1<- x[7] %>% as.numeric
    allele_num1<- x[9] %>% as.numeric
    if(allele_num1 == 1){
      return(NA)
    }else{
      allele_dat <- gene_dat_annote %>% filter(chrome == chrome1, position == position1)
      allele_dat1 <- allele_dat[,-c(1:3)] %>% t() 
      allele_dat2 <- allele_dat1 %>% as.data.frame %>% mutate(regist_no = rownames(allele_dat1))
      colnames(allele_dat2)[1] <- "allele"
      sub_dat <- sub_methyl %>% left_join(allele_dat2, by = "regist_no") %>% filter(allele != "NN")
      kruskal_smry <- kruskal.test(methyl ~ allele, data = sub_dat) 
      p_v<- kruskal_smry$p.value
      return(p_v)
    }
  }
  
  sub_result_new$anova <- apply(sub_result_new, 1, anova_p)
  sub_result_new$kruskal <- apply(sub_result_new, 1, kruskal_p)
  
  result <- rbind(result, sub_result_new)
  
  print(i)
}


# fwrite(result, "EWAS/MDS_EWAS/report_0612/meqtl_data.csv", row.names = F)

### Box plot --------------------------
library(ggplot2)

result_sig_anova <- result %>% filter(anova <0.05)
result_sig_kruskal <- result %>% filter(kruskal <0.05)

dim(result_sig_anova)
dim(result_sig_kruskal)

EWAS_boxplot<- function(data){
  n <- nrow(data)
  plot_list <- list()
  for(i in 1:n){
    
    # Methylation dat processing
    methyl_n <- data$cg_site[i]
    methyl_chrome <- data$cg_chrome[i]
    methyl_position <- data$cg_position[i]
    
    age <- data$year[i]
    if(age == 2){
      bmiq_dat <- bmiq_data_2yr_sel
    }else{
      bmiq_dat <- bmiq_data_6yr_sel
    }
    bmiq_dat1 <- bmiq_dat %>% filter(cg_site == methyl_n) %>% t()
    bmiq_dat2 <- bmiq_dat1[-1, ] %>% as.data.frame
    bmiq_dat3 <- data.frame(methyl = bmiq_dat2$. %>% as.numeric, regist_no = bmiq_dat2 %>% rownames)
    methyl_lab = paste0("DNA methylation of \n", methyl_n, " (", age, "yr)")
    
    # SNP dat processing 
    SNP_n <- data$SNP[i]
    SNP_chrome <- data$SNP_chrome[i]
    SNP_position <- data$SNP_position[i]
    
    SNP_dat <- gene_dat_annote %>% filter(SNP == SNP_n, chrome == SNP_chrome, position == SNP_position) %>% t()
    SNP_dat1 <- SNP_dat[-c(1:3),] %>% as.data.frame
    SNP_dat2 <- data.frame(SNP = SNP_dat1$., regist_no = rownames(SNP_dat1))
    
    if(SNP_n == "---") SNP_n == "Uncharacterized"
    # SNP_lab = paste0(SNP_n, " (", SNP_chrome, ":", SNP_position, ")")
    SNP_lab = paste0(SNP_n)
    
    
    plot_dat <- bmiq_dat3 %>% left_join(SNP_dat2, by = "regist_no")
    plot_dat_omitNN <- plot_dat %>% filter(SNP != "NN")
    pp <- ggplot(plot_dat_omitNN, aes(y = methyl, x = SNP, fill = SNP))+geom_boxplot()+ xlab(SNP_lab) + ylab(methyl_lab) +
      geom_jitter(shape=16, position=position_jitter(0.05), size = 2) + scale_fill_brewer(palette = "Set1") + 
      theme(legend.position = "none", text = element_text(size = 40), axis.text.x = element_text(size = 45))
    plot_list[[i]] <- pp
  }
  return(plot_list)
}

plot_anova <- EWAS_boxplot(data = result_sig_anova)
plot_kruskal <- EWAS_boxplot(data = result_sig_kruskal)

length(plot_anova)
length(plot_kruskal)


pdf("EWAS/MDS_EWAS/report_0612/anova_boxplot.pdf")
for(i in 1:length(plot_anova)){
  pp <- plot_anova[[i]]
  print(pp)
}
dev.off()


pdf("EWAS/MDS_EWAS/report_0612/kruskal_boxplot.pdf")
for(i in 1:length(plot_kruskal)){
  pp <- plot_kruskal[[i]]
  print(pp)
}
dev.off()

