#In this part, we leveraged the annotated scRNA-seq to generate a reference profile, and then deconvolve TCGA bulk RNA-seq with the reference profile to obtain cellular composition in TME.
#To examine the robustness of CIBERSORTx, we applied two additional deconvolution methods, MuSiC and Bisque to compare the results from CIBERSORTx.

#########1.Run CIBERSORTx#############
#We adpoted to run CIBERSORTx on our local server, to achieve this purpose, you need to obtain the token from CIBERSORTx website (https://cibersortx.stanford.edu/).
systemctl start docker
systemctl status docker
docker run -v /home/data/immunQTL/1-cibersortx:/src/data -v /home/data/immunQTL/1-cibersortx/output:/src/outdir cibersortx/fractions --username XXX --token XXX --single_cell TRUE --refsample single_ref_count.txt --mixture TCGA_mRNA_counts.txt -rmbatchSmode TRUE --perm 500 --verbose TURE --G.min 300 --G.max 500 --q.value 0.01 --fraction 0


#########2.Run Music#############
# install devtools if necessary
install.packages('devtools')
# install the MuSiC package
devtools::install_github('xuranw/MuSiC')
# load
library(Seurat)
library(SeuratObject)
library(MuSiC)
library(SingleCellExperiment)

# Data Preparations
# Bulk data
data_df_bulk <- "~/BLCA_mRNA_counts.txt"
BLCA.exprs <- as.matrix(read.table(data_df_bulk, header=TRUE, sep="\t",
                                   row.names=1,
                                   as.is=TRUE))
minimalSet <- ExpressionSet(assayData=BLCA.exprs)

# single cell dataset
load("~/sce.BLCA_anno.Rdata")
BLCA.sce <- SingleCellExperiment(
  assays = list(counts = sce.BLCA@assays[["RNA"]]@counts ), 
  colData =   sce.BLCA@meta.data
)

# Estimate cell type proportions
Est.prop.BLCA = music_prop(bulk.mtx = BLCA.exprs,
                             sc.sce = BLCA.sce, 
                             clusters = 'Immune_All_Low',
                             samples = 'sample_id',
                             verbose = F)

names(Est.prop.BLCA)
#[1] "Est.prop.weighted" "Est.prop.allgene"  "Weight.gene"       "r.squared.full"    "Var.prop"  
# Est.prop.weighted: data.frame of MuSiC estimated proportions, subjects by cell types;
# Est.prop.allgene: data.frame of NNLS estimated proportions, subjects by cell types;
# Weight.gene: matrix, MuSiC estimated weight for each gene, genes by subjects;
# r.squared.full: vector of R squared from MuSiC estimated proportions for each subject;
# Var.prop: matrix of variance of MuSiC estimates.

# save Estimation of cell type proportions
write.csv(Est.prop.BLCA[["Est.prop.weighted"]], file = "BLCA_prop_MuSic.csv", row.names = TRUE)


#########3.Compare the deconvolution results between CIBERSORTx and Music#############
prop_CIBERSORT <- read.csv("BLCA-prop-CIBERSORTx.csv", header = TRUE, row.names = 1)
prop_MuSic <- read.csv("BLCA_prop_MuSic.csv", header = TRUE, row.names = 1)

all(rownames(prop_MuSic) == rownames(prop_CIBERSORT))

cor_matrix_celltype <- c()
for (i in seq(1,12)){
  cor_matrix_celltype <- c(cor_matrix_celltype, cor(prop_CIBERSORT_filter_DCs[,i], prop_MuSic_filter_DC[,i]))
}
print(cor_matrix_celltype)


#########4.Run Bisque#############
#Data reading
load("~/sce.LUAD.Rdata")
LUAD_bk <- read.table("~/TCGA-LUAD_mRNA_counts_19938.txt", header = TRUE,row.names = 1, sep = "\t")
LUAD_cb <- read.table("~/CIBERSORTx_Results.txt", header = TRUE,row.names = 1, sep = "\t")
Idents(sce.LUAD) <- "Immune_All_Low"
single.cell.expression.set <- SeuratToExpressionSet(sce.LUAD, delimiter='_', position=1, version="v3")
LUAD_bk<-as.matrix(LUAD_bk)
bulk.expression.set <- ExpressionSet(assayData = LUAD_bk)
LUAD_res <- ReferenceBasedDecomposition(bulk.expression.set, single.cell.expression.set,use.overlap = FALSE)
names(LUAD_res)

#[1] "bulk.props"       "sc.props"         "rnorm"            "genes.used"       "transformed.bulk"
#bulk.props:a matrix of cell type proportion estimates with cell types as rows and individuals as columns;
#sc.props: a matrix of cell type proportions estimated directly from counting single-cell data;
#rnorm: Euclidean norm of the residuals for each individual's proportion estimates
#genes.used:vector of genes used in decomposition
#transformed.bulk:contains the transformed bulk expression used for decomposition.



#########3.Compare the deconvolution results between CIBERSORTx and Bisque#############
LUAD_bis <- t(LUAD_res$bulk.props)
LUAD_cor <- diag(cor(LUAD_bis,LUAD_cb,method = "pearson"))
LUAD_mid <- median(LUAD_cor)
LUAD_mean <- rowMeans(LUAD_bis)


