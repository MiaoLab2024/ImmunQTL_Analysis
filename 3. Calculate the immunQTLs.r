#In this part, we assessed the genetic effect on cellular composition with the adjustment of covariates, such as principal components, PEER factors, patients’ age, sex and tumor stage. 

#########1.Run PCA to obtain principal components#############
plink --bfile chr_all_imputed_qc_new2 --keep genotype_ID.txt --maf 0.05 --geno 0.05 --hwe 0.000001 --make-bed --out imputed_immunqtl2
plink --bfile imputed_immunqtl2 --indep-pairwise 50 50 0.25 --out imputed_immunqtl2
plink --bfile imputed_immunqtl2 --extract imputed_immunqtl2.prune.in --pca 10 --out imputed_immunqtl2

#########2.Calculate PEER factors#############
library(peer)
simple_unsupervised_demo <- function(){
  y = read.csv("residuals_for_peer.csv",header=T)
  K = 10
  Nmax_iterations = 1000
  model = PEER()
  PEER_setNk(model, K) 
  PEER_setPhenoMean(model, as.matrix(y)) 
  PEER_setPriorAlpha(model,0.001,0.1);
  PEER_setPriorEps(model,0.1,10.);
  PEER_setNmax_iterations(model,Nmax_iterations)
  PEER_update(model)
  X = PEER_getX(model)
  W = PEER_getW(model)
  Alpha = PEER_getAlpha(model)
  Yc = PEER_getResiduals(model)
  results<-data.frame(factor=X )
  write.csv(results," peer_factor.csv")
  pdf("r_demo.pdf",width=8,height=8)
  plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance", main="")
  dev.off()
}
simple_unsupervised_demo()

#########3.Calculate immunQTL#############
library(MatrixEQTL)
useModel = modelLINEAR
SNP_file_name = paste("imputed_immunqtl2.traw")
expression_file_name = paste("CIBERSORTx_Results_filterd.csv")
covariates_file_name = paste("covariates.csv")
output_file_name = tempfile()
pvOutputThreshold = 0.05
errorCovariance = numeric()
snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 2000
snps$LoadFile(SNP_file_name)
gene = SlicedData$new()
gene$fileDelimiter = ","
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 5000
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new()
cvrt$fileDelimiter = ","
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$LoadFile(covariates_file_name)
me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)



