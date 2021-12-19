###################################################
#  먼저 필요한 패키지를 모두 다운받으셔야 합니다.
# install.packages("data.table")
# install.packages("dplyr")
# install.packages("diptest")
# install.packages("stringr")
# install.packages("foreach")
# install.packages("bit64")
# install.packages("doParallel")
###################################################

library(data.table)
library(dplyr)
# library(diptest)
library(stringr)
library(foreach)
# library(doParallel)

# setwd("//DESKTOP-IADPKG9/work/EWAS_EDC")
setwd("/media/kevin/SNU_HD/work/EWAS_EDC")

#### methylation 데이터를 불러오는 코드입니다.  ------------
####데이터가 2세, 6세 각각 bmiq normalization을 한 것과 안한것 이렇게 총 4개가 있습니다. 

bmiq_data_2yr <- fread("cg_data/2yr/bmiq.norm.data.txt", header = T)
bmiq_data_6yr <- fread("cg_data/6yr/bmiq.norm.data.txt", header = T)
# raw_data_2yr <- fread("F:/작업/윤정쌤_methyl/2yr/data.txt", header = T)
# raw_data_6yr <- fread("F:/작업/윤정쌤_methyl/6yr/data.txt", header = T)

### the number of samples ---------------
n_6yr <- 80
n_2yr <- 60

#### 그런데 bmiq_data가 아니라 raw 데이터 같은 경우에는 변수가 다른 값들이 추가 되어있습니다. 
#### 그래서 데이터를 정제해주는 함수를 만들었습니다. raw data를 넣으시면 bmiq 데이터처럼 변수이름이랑 methyl 변수만 추려줍니다.
#### bmiq data를 넣으면 그대로 나오기 때문에 넣으셔도 됩니다. 단 정확한 샘플 숫자를 기재해야 합니다.
extract_methyl <- function(data, n){
  if(ncol(data)>(n+1)){
    ss <- str_split(colnames(data), "[.]")
    col1 <- sapply(ss, function(x) ifelse(x[2] == "AVG_Beta", 1, 0))
    col_candi <- c(1,which(col1 == 1))
    data2 <- data[,..col_candi]
    colnames(data2) <- colnames(data2) %>% str_split("[.]") %>% sapply(function(x) x[1])
  }else{
    data2 <- data
  }
  return(data2)
}
# raw_data_2yr_m <- extract_methyl(raw_data_2yr, n_2yr)
bmiq_data_2yr_m <- extract_methyl(bmiq_data_2yr, n_2yr)
bmiq_data_6yr_m <- extract_methyl(bmiq_data_6yr, n_6yr)

#### edc 데이터를 불러오는 코드입니다. --------------
edc_raw <- fread("edc_data/edc_raw.csv", header = T, fill = T)

### cell count 데이터를 불러들어옵니다. -----------
cell_2yr <- fread("estimate_cell_count/result/cell_2yr_epic.csv", header = T, fill = T)
cell_6yr <- fread("estimate_cell_count/result/cell_6yr_450k.csv", header = T, fill = T)
colnames(cell_2yr)[2] <- "regist_no"
colnames(cell_6yr)[2] <- "regist_no"
## 필요한 cell count 변수만 추려냅니다. ------------
cell_2yr_f <- cell_2yr %>% select(regist_no, CD8T, CD4T, NK, Bcell, Mono, Neu)
cell_6yr_f <- cell_6yr %>% select(regist_no, CD8T, CD4T, NK, Bcell, Mono, Neu)

## 분석에 필요한 변수를 정의해줍니다. ---------------------------
## 이 작업은 값들 중에 숫자가 아닌 이상한 값이 있어서 진행하는 것입니다. 
## 만약 다 숫자 값만 존재한다면 나이변수 만들기 전까지 skip 하셔도 됩니다.
mother_var <- "BPA_M_creatinine MEHHP_M_creatinine MEOHP_M_creatinine MnBP_M_creatinine MECCP_M_cre MBzP_M_cre_2 M3_PBA_M_creatinine Pb_M Hg_M Cd_M Mn_M"
mother_var1 <- unlist(str_split(mother_var, " "))

child_var <- "BPA_creatinine BPS_cre BPF_cre MEHHP_creatinine MEOHP_creatinine MnBP_creatinine MECCP_creatinine MBzP_creatinine_2 b3_PBA_creatinine Pb_C Hg_C Cd_C Mn_C"
child_var1 <- unlist(str_split(child_var, " "))

child_outcome <- "b_wisc_1 scqtotal ars_total1 ars_total2 ars_total"
child_outcome1 <- unlist(str_split(child_outcome, " "))

select_coln <- which(colnames(edc_raw) %in% c(mother_var1, child_var1, child_outcome1))
# select_coln %>% length()
## 위에 있는 변수들은  "검체없음" 같은 요상한 값들이 많아서 다 공백으로 처리하였습니다.
# warnings 은 무시하셔도 됩니다. 그냥 다 숫자 아닌것들은 결측처리 했다고 나오는 경고문입니다.
edc_numeric <- mutate_all(edc_raw[,..select_coln], function(x) as.numeric(as.factor(as.character(x))))
str(edc_numeric)

## 바꾼 수치형 변수를 기존의 key 변수와 합치는 작업입니다.----------
edc_final <- cbind(edc_raw[,-..select_coln], edc_numeric)
# which(colnames(edc_final) %in% c(mother_var1, child_var1, child_outcome1)) %>% length()

## 변수 잘 들어갔는지 확인하는 코드입니다.-------
dim(edc_final);dim(edc_raw)



# which(colnames(edc_raw) == "Smoothed_creatitine")
# which(colnames(edc_raw) == "regist_no")

## create old variable ---------------
## 위의 단계 (숫자 값으로 바꾸는 과정) skip 한경우 밑의 코드를 돌려주세요
# edc_final <- edc_raw

## 나이를 계산해봤습니다.
edc_final$old <- apply(edc_final$age %>% as.data.frame() ,1,
                       function(x)ifelse(x<30, 2, ifelse(x<57, 4, ifelse(x<80, 6, 8))))

## 나이별로 샘플 수를 계산한 것입니다.
edc_2yr <- edc_final %>% subset(old == 2);edc_2yr %>% nrow()
edc_4yr <- edc_final %>% subset(old == 4);edc_4yr %>% nrow()
edc_6yr <- edc_final %>% subset(old == 6);edc_6yr %>% nrow()
edc_8yr <- edc_final %>% subset(old == 8);edc_8yr %>% nrow()

# edc_m_6yr <- edc_m %>% subset(old == 6)
## edc 데이터와 cell_count 데이터를 합칩니다. 이때 샘플수도 60, 80으로 맞춰집니다.
edc_2yr_m <-merge(cell_2yr_f, edc_2yr, by = "regist_no", all.x = T)
edc_2yr_m %>% nrow()
edc_6yr_m <-merge(cell_6yr_f, edc_6yr, by = "regist_no", all.x = T)
edc_6yr_m %>% nrow()


EWAS_linear <- function(methyl_data, edc_data, outcome, predictor, covariate = c(),
                        n = nrow(methyl_data)){
  st_t <- Sys.time()
  # n <- nrow(methyl_data)
  # n <- 10^5
  cat("Start:", paste(Sys.time(), "\n"))
  # candi <- data.frame(estimate = numeric(n),
  #                     p_value = numeric(n),
  #                     rownumber = integer(n))
  nc <- ncol(methyl_data)
  if(outcome == "methyl"){
    edc_data_select <- edc_data %>% select(regist_no, predictor, covariate, CD8T, CD4T, NK, Bcell, Mono, Neu)
  }else if(predictor == "methyl"){
    edc_data_select <- edc_data %>% select(regist_no, outcome, covariate, CD8T, CD4T, NK, Bcell, Mono, Neu)
  }

  fit_linear <- function(i, methyl_data, outcome, predictor, covariate){
    sub <- t(methyl_data[i,2:nc])
    regist_no <- as.numeric(rownames(sub))
    sub <- cbind(sub, regist_no)
    rownames(sub) <- NULL
    colnames(sub)[1] <- "methyl"
    sub_d <- merge(edc_data_select, sub, by ="regist_no", all.y = T)
    # sub_d1 <- select(sub_d, outcome, predictor, covariate)
    sub_d1 <- sub_d %>% select(outcome, predictor, covariate)
    formula1 <- as.formula(paste0(outcome," ~ ."))
    tryCatch({fit <- summary(lm(formula1, data = sub_d1))
    })
    # fit <- summary(lm(formula1, data = sub_d1))
    result <- cbind(i, fit$coefficients[2,4], fit$coefficient[1,4])
    return(result)
  }
  # progress <- function(n) 
  # opts <- list(progress = progress)
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)

  progcombine<-function(){
    Sys.sleep(0.1)
    pb <- txtProgressBar(min=1, max=n-1,style=3)
    count <- 0
    function(...) {
      count <<- count + length(list(...))-1
      setTxtProgressBar(pb,count)
      flush.console()
      rbind(...)
    }
  }
  candi1 <- foreach(i = icount(n), .combine = progcombine(), .packages = c('dplyr')) %dopar%{
    # setTxtProgressBar(k, j)
    # Sys.sleep(0.1)
    s<- fit_linear(i, methyl_data, outcome, predictor, covariate)
    return(s)
    # setTxtProgressBar(pb, i)
    # regist_no <- sub %>% rownames() %>% as.numeric()
    
    
    # if(na.omit(sub_d1) %>% nrow() > 2 ){ 

      # fit <- lm(formula1, data = sub_d1) %>% summary()
      
      # }
      # error=function(e){cat("ERROR in", i,"th row :",conditionMessage(e), "\n")})
      # if(i %% 10^3 == 0){
      # est_t <- ((i/n-1)*(st_t %>% as.numeric())+(Sys.time() %>% as.numeric()))*(n/i)
      # cat("processing ", i, "of", n," : ",
      #      format(round(i/n*100, digits = 1), nsmall = 1),
      #     "%.", "| estimated end time: ", paste(as.POSIXct(est_t, origin = "1970-01-01")), "\n")
      # }
    # }
          }
  # close(pb)
  stopCluster(cl)
  # candi$p_value <- candi$p_value %>% as.character() %>% as.numeric()
  total_cg <- t(methyl_data[(1:n),1]) %>% as.character()
  candi <- candi1 %>% as.data.frame
  colnames(candi) <- c("rownumber", "p_value", "estimate")
  candi$cg_site <- total_cg[candi$rownumber]
  select_ind <- which(candi$estimate <= 0)
  if(select_ind %>% length == 0){
    final_candi <- candi
  }else{
    final_candi <- candi[-select_ind, ]
  }
  # cg_select <- total_cg[candi$rownumber]
  # final_candi$cg_site <- cg_select[pv_index]
  rownames(final_candi) <- NULL
  cat("Total processing time: ",
      paste(difftime(Sys.time(), st_t, units = 'mins')) %>% as.numeric() %>%  round(digits = 2),
      "mins", "\n")
  return(final_candi)
  
}


### 제가 함수를 만들었습니다. 최대한 필요한 변수를 넣을 수 있게 만들었습니다.
# 1. methyl_data 는 methylation 데이터를 넣으시면 됩니다. 
# 2. edc_data는 정제된 edc 데이터를 넣으면 됩니다.
# 3. outcome(반응변수) 과 predictor(독립변수) 중에 methylation을 넣고 싶으신 곳에 "methyl" 이라고 넣으시면 됩니다.
# 3. 나머지 outcome 이나 predictor에는 원하시는 edc변수를 넣으시면 됩니다. 변수는 edc 데이터에 있는 변수명과 일치해야합니다.
# 4. method는 pvalue를 보정하는 방법입니다. 코드창에 ?p.adjust라고 치면 나오는 방법들을 모두 넣을 수 있습니다.
# 5. 조절변수는 covariate로 넣으시면 되고 제가 해봤는데 baby_gen(성별)은 가끔 na 값을 지우면 male 이나 female 만 남는 경우가 있어서 오류가 납니다.
# 그래서 저는 키와 몸무게만 했으며, covariate에 들어갈 변수 역시 edc 데이터에 정확한 변수명으로 존재해야합니다.
# 그리고 covariate에 넣으실때 c() 안에 "" 표시하셔서 넣으셔야 합니다.
# 6. 결측값을 지웠을 때 관측값이 2개 이하인 경우 skip 합니다.

library(doParallel) 
library(doSNOW)
library(tcltk)

mother_var1;child_var1
for(j in 1:2){
  edc_data = c("edc_2yr_m", "edc_6yr_m")
for(i in 1:length(child_var1)){
test <- EWAS_linear(methyl_data = bmiq_data_6yr_m, edc_data = edc_data[j] %>% get,
                    outcome = "methyl", predictor = child_var1[i])

file_n <- paste0("/media/kevin/SNU_HD/work/EWAS_EDC/result/methyl_6yr_",child_var1[i],"_",edc_data[j],"_crude.csv")
write.csv(test, file_n)
}
}


test_6cov <- EWAS_linear(methyl_data = raw_data_2yr_m, edc_data = edc_2yr_m,
                    outcome = "methyl", predictor = "MEHHP_creatinine_2", covariate = c("bmi", "baby_gen",
                                                                                        "fa_hei",	"fa_wei", "mo_hei",	"mo_wei", "MEHHP_M_creatinine_2"
                    ))
# test_2cov[which(test_2cov$rownumber == 1366),]

write.csv(test_6cov, "F:/작업/윤정쌤_methyl/EWAS/result/single_MEHHP_creatinine_2_6cov.csv")


## test에는 cg_site, estimate(베타값), p_value, 그리고 rownumber가 있습니다. 
## rownumber는 해당 cg_site의 행 위치입니다.



