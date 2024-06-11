library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(clustree)
library(cowplot)

#########1.Create the Seurat object##########
dir='/home/data/0-data' 
samples=list.files( dir )
samples 
#Load data
sceList<-list()
for (i in 1:length(samples)){
  counts <- Read10X_h5(samples[i])
  sceList[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200,project = samples[i])
}

#########2.Quality control#############
#Calculate percent ribosomal and mitochondrial genes 
for (i in 1:length(sceList)){
  sceList[[i]][["percent.mt"]]<-PercentageFeatureSet(sceList[[i]], pattern = "^MT-")
  sceList[[i]][["percent.rb"]]<-PercentageFeatureSet(sceList[[i]], pattern = "^RP[SL]")
  sceList[[i]][["percent.HB"]]<-PercentageFeatureSet(sceList[[i]], pattern ="^HB[^(P)]")
}  

#Filter cells based on features, count and percent MT 
library(scater)
for (i in 1:length(sceList)){
  sceList[[i]]@meta.data$nCount_RNA_outlier_2mad<-isOutlier(log(sceList[[i]]@meta.data$nCount_RNA),log = F,type = "lower",nmads = 2)
  sceList[[i]]@meta.data$nFeature_RNA_outlier_2mad<-isOutlier(log(sceList[[i]]@meta.data$nFeature_RNA),log = F,type = "lower",nmads = 2)
  sceList[[i]]@meta.data$percent_mt_outlier_2mad<-isOutlier(log(sceList[[i]]@meta.data$percent.mt),log = F,type = "lower",nmads = 2)
}  

for (i in 1:length(sceList)){
  sceList[[i]]<-subset(sceList[[i]],subset=nCount_RNA_outlier_2mad == "FALSE" & 
                         nFeature_RNA_outlier_2mad == 'FALSE' &
                         percent_mt_outlier_2mad == "FALSE") 
  sceList[[i]]<- NormalizeData(sceList[[i]],assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000)
  
}

for (i in 1:length(sceList)) {
  sceList[[i]] <- NormalizeData(sceList[[i]],assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
  sceList[[i]] <- FindVariableFeatures(sceList[[i]], selection.method = "vst",nfeatures = 2000)
  sceList[[i]] <- ScaleData(sceList[[i]], vars.to.regress = c("nFeature_RNA", "percent.mt"))
  sceList[[i]] <- RunPCA(sceList[[i]], pc.genes = VariableFeatures(sceList[[i]]))
  sceList[[i]] = RunTSNE(sceList[[i]], dims = 1:15)
  sceList[[i]] = RunUMAP(sceList[[i]], dims = 1:15)
  sceList[[i]] <- FindNeighbors(sceList[[i]], dims = 1:15)
  sceList[[i]] <- FindClusters(sceList[[i]], resolution = 0.8)
}

#remove doulet
library(DoubletFinder)
for (i in 1:length(sceList)) {
  sweep.res.list <- paramSweep_v3(sceList[[i]], PCs = 1:10, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats) 
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() 
  DoubletRate1 = ncol(sceList[[i]])*8*1e-6 
  homotypic.prop <- modelHomotypic(sceList[[i]]$seurat_clusters) 
  nExp_poi <- round(DoubletRate1*ncol(sceList[[i]])) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  sceList[[i]] <- doubletFinder_v3(sceList[[i]], PCs = 1:10, pN = 0.25, pK = pK_bcmvn, 
                                   nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  colnames(sceList[[i]]@meta.data)[ncol(sceList[[i]]@meta.data)]="DoubletFinder"
  table(sceList[[i]]@meta.data$DoubletFinder)
  singlet_cells <- row.names(sceList[[i]]@meta.data)[which(sceList[[i]]@meta.data$DoubletFinder=='Singlet')]
  sceList[[i]]=subset(sceList[[i]],cells=singlet_cells)
}
scRNA1<-merge(sceList[[1]],sceList[2:length(sceList)],add.cell.ids = samples) 


#########3.Normalize and scale the data#############
scRNA1 <- NormalizeData(scRNA1,assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst", nfeatures = 2000)
scRNA1 <- ScaleData(scRNA1, vars.to.regress = c("nCount_RNA", "percent.mt"))
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1))
ElbowPlot(scRNA1, ndims = 50,reduction = "pca")
pc.num=15
library(harmony)
scRNA1 <- RunHarmony(scRNA1, group.by.vars = "orig.ident")
names(scRNA1@reductions)
harmony_embeddings <- Embeddings(scRNA1, 'harmony')
scRNA1 <- RunTSNE(scRNA1, reduction = "harmony", dims = 1:15)
scRNA1 <- RunUMAP(scRNA1, reduction = "harmony", dims = 1:15)
scRNA1=FindNeighbors(scRNA1, reduction = "harmony", dims = 1:15)
scRNA1=FindClusters(scRNA1,resolution = 0.8)
p1=DimPlot(scRNA1,reduction = "umap",group.by = "RNA_snn_res.0.8",label=T)
p2=DimPlot(scRNA1,reduction = "umap",group.by = "orig.ident",label=T)
plota=p1+p2

#########4.General annotations#############
celltype_marker=c(
  "EPCAM","KRT18","KRT17","KRT14",#Epithelial cells
  "THY1","COL1A1","COL1A2","FGF7",#Fibroblasts
  "PECAM1","CDH5","VWF","ENG",#Endothelial cells
  "PTPRC","CD3E","CD14","CD79A","MS4A2"#Immune cells
)
DotPlot(scRNA1, features = celltype_marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

#Annotate Immune vs Nonimmune clusters
#At this point we dont care for a more detailed annotation as we will annotate immune cells separately later
celltype=data.frame(ClusterID=0:17,
                    celltype=0:17) 
celltype[celltype$ClusterID %in% c(5,11),2]='Epithelial cells' 
celltype[celltype$ClusterID %in% c(10,16),2]='Fibroblasts' 
celltype[celltype$ClusterID %in% c(8),2]='Endothelial cells' 
celltype[celltype$ClusterID %in% c(0:4,6,7,9,12:15,17),2]='Immune cells' 
for(i in 1:nrow(celltype)){
  scRNA1@meta.data[which(scRNA1@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(scRNA1@meta.data$celltype)
p3=DimPlot(scRNA1, reduction = "umap", group.by = "celltype",label = T) 
p4=DimPlot(scRNA1, reduction = "umap", group.by = "RNA_snn_res.0.8",label = T) 
plotb=p3+p4
ggsave('Allcell_umap_by_celltype.png',device = "png")

#Generate immune and nonimmune cell lists
cells.use.immune <- row.names(scRNA1@meta.data)[which(scRNA1@meta.data$celltype=='Immune cells')]
scRNA1.CD45neg<-subset(scRNA1, cells=cells.use.immune,invert=T)
scRNA1.CD45pos<-subset(scRNA1, cells=cells.use.immune)

#########5.Annotate immune cells with CellTypist#############
library(harmony)
scRNA1.CD45pos <- NormalizeData(scRNA1.CD45pos,assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000)
scRNA1.CD45pos <- FindVariableFeatures(scRNA1.CD45pos, selection.method = "vst", nfeatures = 2000)
scRNA1.CD45pos <- ScaleData(scRNA1.CD45pos, vars.to.regress = c("nFeature_RNA", "percent.mt"))
scRNA1.CD45pos <- RunPCA(scRNA1.CD45pos, features = VariableFeatures(scRNA1.CD45pos))
ElbowPlot(scRNA1.CD45pos, ndims = 50,reduction = "pca")

scRNA1.CD45pos <- RunHarmony(scRNA1.CD45pos, group.by.vars = "orig.ident")
names(scRNA1.CD45pos@reductions)
harmony_embeddings <- Embeddings(scRNA1.CD45pos, 'harmony')
#harmony_embeddings[1:5, 1:5]

scRNA1.CD45pos <- RunTSNE(scRNA1.CD45pos, reduction = "harmony", dims = 1:20)
scRNA1.CD45pos <- RunUMAP(scRNA1.CD45pos, reduction = "harmony", dims = 1:20)
scRNA1.CD45pos=FindNeighbors(scRNA1.CD45pos, reduction = "harmony", dims = 1:20)
scRNA1.CD45pos=FindClusters(scRNA1.CD45pos,resolution = 0.8)
p5=DimPlot(scRNA1.CD45pos,reduction = "umap",group.by = "RNA_snn_res.0.8",label=T)
p6=DimPlot(scRNA1.CD45pos,reduction = "umap",group.by = "orig.ident",label=T)
plotc=p5+p6
ggsave(plot=plotc,"Immune_UMAP_after_harmony.png")

#Annotate with marker genes
celltype_marker2=c("CD3D",'CD3E','CD2',"CD4","CD8A",#T cell
                   "IFNG","NKG7",#Teff 
                   "GZMB","PRF1",#Tcyto
                   "BCL6","CXCR5",#Tfh
                   "IFNG","GZMB",#Th1
                   'CD79A','MZB1','MS4A1','CD79B',#B cell
                   'FOXP3',"IL32",'TNFRSF18','TNFRSF4',#Treg
                   "LYZ",'SEPP1','C1QA','APOE','CD14','CD68','RNASE1',#Macrophage
                   'TPSAB1','TPSB2','CPA3','HPGDS',#Mast
                   'HLA-DRA','HLA-DPB1','CST3','HLA-DPA1',#mDC
                   'PTGDS','SOX4','GZMB','IRF7',#pDC
                   'IGHA1','IGHG1',"IGHG2",#Plasma
                   "NKG7",'KLRF1','KLRD1'#NK
)
DotPlot(scRNA1.CD45pos, features = celltype_marker2)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave(filename = "Immune_marker.png",device = "png",width = 44,height = 33,units = "cm")

celltype=data.frame(ClusterID=0:20,
                    celltype=0:20) 
celltype[celltype$ClusterID %in% c(0,11,12,15,16,18),2]='Tcell' 
celltype[celltype$ClusterID %in% c(5),2]='Treg' 
celltype[celltype$ClusterID %in% c(2,7,13),2]='Macrophage' 
celltype[celltype$ClusterID %in% c(4,6,19,20),2]='Bcell' 
celltype[celltype$ClusterID %in% c(1,3,10),2]='Teff' 
celltype[celltype$ClusterID %in% c(8),2]='Tcyto' 
celltype[celltype$ClusterID %in% c(17),2]='Th1' 
celltype[celltype$ClusterID %in% c(9),2]='NK' 
celltype[celltype$ClusterID %in% c(14),2]='Mast cell' 

scRNA1.CD45pos@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA1.CD45pos@meta.data[which(scRNA1.CD45pos@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(scRNA1.CD45pos@meta.data$celltype)
p7=DimPlot(scRNA1.CD45pos, reduction = "umap", group.by = "celltype",label = T) 
p8=DimPlot(scRNA1.CD45pos, reduction = "umap", group.by = "RNA_snn_res.0.8",label = T) 
plotd=p7+p8
ggsave('Immune_harmony_umap_by_celltype.png',device = "png")
save(scRNA1.CD45pos,file = "scRNA1.CD45pos.harmony.Rdata")

#Covert Seurat to H5ad
library(Seurat)
library(SeuratDisk)
library(SeuratData)
load("scRNA1.CD45pos.harmony.Rdata")
SaveH5Seurat(scRNA1.CD45pos,"scRNA1.CD45pos.harmony.h5Seurat")
Convert("scRNA1.CD45pos.harmony.h5Seurat",dest="h5ad")

#Annotate with CellTypist
library(reticulate)
reticulate::use_condaenv(condaenv = "scanpy", required = TRUE)
scanpy = import("scanpy")
celltypist = import("celltypist")
pandas <- import("pandas")
numpy = import("numpy")
numba=("numba")
adata1=scanpy$read_h5ad("scRNA1.CD45pos.harmony.h5ad")

model = celltypist$models$Model$load(model = 'Immune_All_Low.pkl')
predictions = celltypist$annotate(adata1, model = 'Immune_All_Low.pkl',majority_voting = T)
predictions$predicted_labels %>% head()
labels<-predictions$predicted_labels
scRNA1.CD45pos  = AddMetaData(scRNA1.CD45pos, predictions$predicted_labels$majority_voting, col.name ="Immune_All_Low")

load("scRNA1.CD45pos.harmony_celltypist.Rdata")
p9 = DimPlot(scRNA1.CD45pos,group.by = "Immune_All_Low", reduction = "umap", label = TRUE,pt.size = 0.5,label.box = T) 
p10 = DimPlot(scRNA1.CD45pos,group.by = "RNA_snn_res.0.8", reduction = "umap", label = TRUE,pt.size = 0.5,label.box = T) 
plote=p9+p10

#Merge CD45neg and CD45pos
scRNA1.CD45neg@meta.data$Immune_All_Low<-scRNA1.CD45neg@meta.data$celltype
sce <- merge(scRNA1.CD45neg, y=scRNA1.CD45pos)

#Remove ribosomal and mitochondrial genes
mt.genes <- rownames(sce)[grep("^MT-",rownames(sce))]
rb.genes <- rownames(sce)[grep("^RP[SL]",rownames(sce))]
ercc.gene<-rownames(sce)[grep("^ERCC-",rownames(sce))] 
del_gene<-c(mt.genes,rb.genes,ercc.gene)
selected_f2<-setdiff(rownames(sce),del_gene)
sce <- subset(sce, features = selected_f2)

single_ref=as.data.frame(GetAssayData(sce, slot='counts')) 
colnames(single_ref)<-sce@meta.data$Immune_All_Low
write.table(single_ref,"single_ref_count.txt",sep = "\t",quote = F)

