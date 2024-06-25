library(ComplexHeatmap)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tidyverse)
library(fgsea)
library(ggplot2)
library(ImmuLncRNA)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)

#########1.Do normalization#############
#Load custom functions
Adjust_mrna <- function(mrna){
  mRow30 <- which(apply(mrna,1,function(v){return((sum(v==0)/length(v))>=0.3)}))
  mRemv <- mRow30
  if(length(mRemv)==0){
      mRNA_out0 <- mrna
    }else{
      mRNA_out0 <- mrna[-(mRemv),]
    }
    mRNA_exp_inter <- log2(mRNA_out0+0.001)
    return(mRNA_exp_inter)
}

gene_expression<-read.table("gene_fpkm_protein_coding_19894gene.txt",header = T,check.names = F,row.names = 1,sep="\t")
selected_columns <- grep("^T", colnames(gene_expression), value = TRUE)
  gene_expression_tumor<-gene_expression[,selected_columns]
  id<-colnames(gene_expression_tumor)
  id <- gsub("T", "", id)
  id1<-sub("_2$", "", id)
  colnames(gene_expression_tumor)<-id1
  acc_mrna <- Adjust_mrna(gene_expression_tumor)

CIBERSORTx<-read.csv("CIBERSORTx_Results_filter.csv",header = T,check.names = F,row.names = 1)
id2<-colnames(CIBERSORTx)
id3<-intersect(id1,id2)
acc_immune<-CIBERSORTx[,id3]
acc_mrna<-acc_mrna[,id3]

#########2.Estimate tumor purity#############
library(estimate)
# Calculates the infiltration level
filterCommonGenes("purityinput.txt", output.f="exp_deal.gct", id="GeneSymbol")
estimateScore("exp_deal.gct", "estimate_score_all.gct", platform="illumina")
est_score_all <- readLines("estimate_score_all.gct")
est_scores <- unlist(strsplit(est_score_all[grep("ESTIMATEScore",est_score_all)],"\t")) %>% 
  .[3:length(.)] %>% as.numeric()
#Calculates tumour purity
Tumour_purity <- cos(0.6049872018+0.0001467884*est_scores)
samples_name <- unlist(strsplit(est_score_all[grep("NAME",est_score_all)],"\t")) %>% 
  .[3:length(.)]
names(Tumour_purity) <- samples_name %>% str_replace_all(., "\\.", "-")

#########3.Calculate partial correlation coefficient#############
fun_mtx_pcr <- function(x,y,z){
    r12=cor(t(x),t(y))
    r13=cor(t(x),z)
    r23=cor(z,t(y))
    r123=r13%*%r23
    rup=r12-r123
    rd1=sqrt(1-r13*r13)
    rd2=sqrt(1-r23*r23)
    rd=rd1%*%rd2
    rm(r12, r13, r23, r123, rd1, rd2)
    gc()
    rrr=rup/rd
    return(rrr)
}

# Extrate samples with both gene expression and immune cell fraction data
  sam <- intersect(colnames(acc_mrna),colnames(acc_immune)) 
  acc_immune<-acc_immune[,sam]
  acc_mrna<-acc_mrna[,sam]
  Tumour_purity<-Tumour_purity[sam]
  
  all(colnames(acc_immune) == colnames(acc_mrna))
  all(colnames(acc_mrna) == names(Tumour_purity))
  n = ncol(acc_immune) #number of samples
  gn = 1 #just one factor need to be corrected: Tumour_purity
  pcor <- fun_mtx_pcr(acc_immune,acc_mrna, Tumour_purity)
  #calculate statistic and p.value
  statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
  p.value <- 2*pnorm(-abs(statistic))
  rownames(pcor) <- rownames(acc_immune) 
  rownames(p.value) <- rownames(acc_immune)
  colnames(pcor) <- rownames(acc_mrna) 
  colnames(p.value) <- rownames(acc_mrna)
  #calculate RS: sort by pvalue and give pvalue a direction based on positive and negative correlation
  RS <- -log10(p.value)*sign(pcor)
  data_melt<-melt (RS)
  names(data_melt) = c('cell_type', 'symbol', 'Rankscore')
  data_p<-melt (p.value)
  names(data_p) = c('cell_type', 'symbol', 'P-value')
  data_pcor<-melt (pcor)
  names(data_pcor) = c('cell_type', 'symbol', 'cor')
  data<-left_join(data_melt,data_p,by=c("cell_type"="cell_type","symbol"="symbol"))
  data<-left_join(data,data_pcor,by=c("cell_type"="cell_type","symbol"="symbol"))
  data<-na.omit(data)
  data$FDR<-p.adjust(data$`P-value`, method="fdr")
  write.csv(data,paste0("gene_immunqtl_cell_correlation_pcc.csv"),quote = F,row.names = F)

