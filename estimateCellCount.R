# 2019.7.31. Yoon-Jung Choi 
# Estimation of cell type fraction 


# Illumina EPIC data on immunomagnetic sorted peripheral adult blood cells
# https://bioconductor.org/packages/release/data/experiment/html/FlowSorted.Blood.EPIC.html
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("FlowSorted.Blood.EPIC")


# To view documentation for the version of this package installed in your system, start R and enter:
browseVignettes("FlowSorted.Blood.EPIC")


# code source:
# http://127.0.0.1:28978/library/FlowSorted.Blood.EPIC/doc/FlowSorted.Blood.EPIC.Rmd

# Step 1: Load the reference library to extract the artificial mixtures  
library(ExperimentHub)  

hub <- ExperimentHub()  

query(hub, "FlowSorted.Blood.EPIC")  

FlowSorted.Blood.EPIC <- hub[["EH1136"]]  

FlowSorted.Blood.EPIC 
dim(FlowSorted.Blood.EPIC)


# Step 2 separate the reference from the testing dataset  

RGsetTargets <- FlowSorted.Blood.EPIC[,
             FlowSorted.Blood.EPIC$CellType == "MIX"]  
             
sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
                            seq_len(dim(RGsetTargets)[2]), sep = "_")  
RGsetTargets  

# Step 3: Deconvolute using the IDOL L-DMR  

head (IDOLOptimizedCpGs)  

# If you need to deconvolute a 450k legacy dataset use 
# IDOLOptimizedCpGs450klegacy instead  

# Do not run with limited RAM the normalization step requires a big amount of 
# memory resources  


#https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylationEPICmanifest.html
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylationEPICmanifest")

library(IlluminaHumanMethylationEPICmanifest)

#위까지가 준비 단계 
#이게 main 분석 

if (memory.limit()>8000){  
 countsEPIC<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood",   
                                processMethod = "preprocessNoob",  
                                probeSelect = "IDOL",  
                                cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),  
                                referencePlatform = "IlluminaHumanMethylationEPIC",  
                                referenceset = NULL,  
                                IDOLOptimizedCpGs =IDOLOptimizedCpGs,   
                                returnAll = FALSE)  
                                
head(countsEPIC$counts)  
}

head(RGsetTargets)
dim(RGsetTargets)
# 1051815 x 12 (12개 mixture(=sample)에서 100만 여개의 행이 있는데 이게 정체가 뭐지? 왜 이렇게 많지??    


head(IlluminaHumanMethylationEPIC)

head(IDOLOptimizedCpGs) 
length(IDOLOptimizedCpGs)
#이거는 450개 cg만 있는데... 

############################################################################################
#idat file 읽어오기 
#illuminaio 설치하기 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("illuminaio")
library(illuminaio)

#IlluminaDataTestFiles 설치하기 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IlluminaDataTestFiles")



d1=read.idat("D:\\2018_EDC\\rawdata\\methylation_data\\2세_60명\\(1)Raw.data\\idat\\201364910042_R01C01_Grn.idat",  
		bgxfile, dateinfo=FALSE, annotation="Symbol", tolerance=0L, verbose=TRUE)
# 안 됨. 

if(require(IlluminaDataTestFiles)) {
  idatFile <- system.file("D:\\2018_EDC\\rawdata\\methylation_data\\2세_60명\\(1)Raw.data\\idat", "idat", "201364910042_R01C01_Grn.idat",
                          package = "IlluminaDataTestFiles")
  idat <- readIDAT(idatFile)
  names(idat)
  idat$Quants[1:5,]
}
# 당연히 안 됨. 

###################################################################################
# The minfi User's Guide
# https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html
# Citing the minfi package 위 사이트 참고 


library(minfi)
library(minfiData)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("minfiData")

# 다음은 minfi 패키지 안에 들어있는 예시다. 
MsetEx <- preprocessRaw(RGsetEx)
GMsetEx <- mapToGenome(MsetEx)


baseDir <- system.file("extdata", package = "minfiData")
list.files(baseDir)
list.files(file.path(baseDir, "5723646052"))

targets <- read.metharray.sheet(baseDir)
targets

#"D:\2018_EDC\rawdata\methylation_data\2세_60명\(1)Raw.data\data.txt"를 SAS에서 불어와서 csv 파일로 변환한 다음 아래와 같이 열었음. 
rgset=read.csv("D:\\2018_EDC\\output\\rgset.csv")
dim(rgset)
# 866297 x 241 
# -> of no use

#  object is of class 'data.frame', but needs to be of class 'RGChannelSet' 'RGChannelSetExtended' or 'MethylSet' to use this function


# RGChannelSet 으로 만들기 
path <- "idats"
list.files(path)
Basename <- file.path("D:\\2018_EDC\\rawdata\\methylation_data\\2세_60명\\(1)Raw.data\\idat\\201364910042_R01C01_Grn.idat", targets$Basename)

read.metharray(Basename, extended=FALSE, verbose = FALSE, force=FALSE)