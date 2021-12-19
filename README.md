# SNU_methylation_coding

* This is the first trial to analyze the methylation dataset.
* I used the beta value of the methylation data
* batcheffect.R is for adjusting the batcheffect of the methylation data. 
* annotaiton.R is for the annotation of the CpG sites
* estimateCellCount.R is for adjusting the cell type of the blood samples. There are several reference to estimate it.
* meqtl.R is for the analyzing mQTL.
* EWAS_linear.R is the linear regression by setting confouonders as above. There are too many linear models to run, so I used parallel computing by foreach function.
