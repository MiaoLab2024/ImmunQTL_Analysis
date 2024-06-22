library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)

##Extract expression matrix and cell classification
data.input <- GetAssayData(sce.BRCA, assay = "RNA", slot = "data")
identity <- subset(sce.BRCA@meta.data, select = "Immune_All_Low_new")


#############################1.cell-cell############################
#Crerat cellchat object
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "Immune_All_Low_new")

####CellChatDB.human and CellChatDB.mouse are optional
dir.create("./cell-cell")
setwd("./cell-cell/")
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
cellchat@DB <- CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat)
#future::plan("multiprocess", workers = 24)
# Identify overexpressed genes
cellchat <- identifyOverExpressedGenes(cellchat)
#Identify ligand-receptor pairs
cellchat <- identifyOverExpressedInteractions(cellchat)
#Map ligands and receptors into the PPI network
cellchat <- projectData(cellchat, PPI.human)

## Interaction inference
## 1、calculate the communication probability to infer the communication network of cell interaction
cellchat <- computeCommunProb(cellchat, raw.use = F,population.size = TRUE) 
###Filter out small cell populations
cellchat <- filterCommunication(cellchat, min.cells = 3)

#Extract inferred cell interaction communication network data matrix, which includes all inferred ligand/receptor level cell-cell communication.
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"net_lr_cell_cell.csv")
levels(cellchat@idents)

##2、Infer intercellular communication at the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) #Count all ligand receptor pairs on a pathway
##Summarize communication probabilities to calculate aggregated communication networks between cells
cellchat <- aggregateNet(cellchat)

##3、Calculate the aggregation cell interaction communication network
groupSize <- as.numeric(table(cellchat@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE) 
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

##############Save results
pathway.show.all=cellchat@netP$pathways
##Generate bubble plot
levels(cellchat@idents)
dir.create("./bubble")
setwd("./bubble/")
for (i in 1:length(levels(cellchat@idents))) {
  p=netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:length(levels(cellchat@idents))), thresh = 0.05,remove.isolate = F,angle.x = 45)
  ggsave(filename = paste0(gsub("/","_",levels(cellchat@idents)[i]),"_bubble.pdf"),
         p,width=8,height=12,dpi=300)
}

##Individual cells interact with other cells
#a.intensity or probability of interactions
setwd("../")
dir.create("./circle")
setwd("./circle/")
dir.create("./weight")
setwd("./weight/")

mat <- cellchat@net$weight
#par(mfrow = c(1,1), xpd=TRUE)
#par("mar")
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  pdf(file = paste0(gsub("/","_",levels(cellchat@idents)[i]),"_circle_weight.pdf"),width=8,height=12)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,vertex.label.cex = 0.5, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

#b.number of interaction
setwd("../")
dir.create("./count")
setwd("./count/")

mat <- cellchat@net$count
#par(mfrow = c(1,1), xpd=TRUE)
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  pdf(file = paste0(gsub("/","_",levels(cellchat@idents)[i]),"_circle_count.pdf"),width=8,height=12)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,vertex.label.cex = 0.5, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}
rm(cellchat)

#########################2.signal################################
setwd("../..")
getwd()
setwd("../")
dir.create("./Secreted_Signaling")
setwd("./Secreted_Signaling/")

#Crerat cellchat object
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "Immune_All_Low_new")

####CellChatDB.human and CellChatDB.mouse are optional
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat)
#future::plan("multiprocess", workers = 24)
# Identify overexpressed genes
cellchat <- identifyOverExpressedGenes(cellchat)
#Identify ligand-receptor pairs
cellchat <- identifyOverExpressedInteractions(cellchat)
#Map ligands and receptors into the PPI network
cellchat <- projectData(cellchat, PPI.human)

## Interaction inference
## 1、calculate the communication probability to infer the communication network of cell interaction
cellchat <- computeCommunProb(cellchat, raw.use = F,population.size = TRUE) 
###Filter out small cell populations
cellchat <- filterCommunication(cellchat, min.cells = 3)

#Extract inferred cell interaction communication network data matrix, which includes all inferred ligand/receptor level cell-cell communication.
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"net_lr_signal.csv")
levels(cellchat@idents)

##2、Infer intercellular communication at the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) #计算一个通路上的所有配体受体对
##Summarize communication probabilities to calculate aggregated communication networks between cells
cellchat <- aggregateNet(cellchat)

##3、Calculate the aggregation cell interaction communication network
groupSize <- as.numeric(table(cellchat@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE) 
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

##############Save results
pathway.show.all=cellchat@netP$pathways
dir.create("./bubbel")
setwd("./bubbel/")

##Generate bubble plot
levels(cellchat@idents)
for (i in 1:length(levels(cellchat@idents))) {
  p=netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:length(levels(cellchat@idents))), thresh = 0.05,remove.isolate = F,angle.x = 45)
  ggsave(filename = paste0(gsub("/","_",levels(cellchat@idents)[i]),"_bubble.pdf"),
         p,width=8,height=12,dpi=300)
}

##Individual cells interact with other cells
#a.intensity or probability of interactions
setwd("../")
dir.create("./circle")
setwd("./circle/")
dir.create("./weight")
setwd("./weight/")

mat <- cellchat@net$weight
#par(mfrow = c(1,1), xpd=TRUE)
#par("mar")
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  pdf(file = paste0(gsub("/","_",levels(cellchat@idents)[i]),"_circle_weight.pdf"),width=8,height=12)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,vertex.label.cex = 0.5, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

#b.number of interaction
setwd("../")
dir.create("./count")
setwd("./count/")

mat <- cellchat@net$count
#par(mfrow = c(1,1), xpd=TRUE)
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  pdf(file = paste0(gsub("/","_",levels(cellchat@idents)[i]),"_circle_count.pdf"),width=8,height=12)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,vertex.label.cex = 0.5, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

#############################3.ECM############################
setwd("../..")
setwd("../")
dir.create("./ECM")
setwd("./ECM/")
rm(cellchat)
#Crerat cellchat object
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "Immune_All_Low_new")

####CellChatDB.human and CellChatDB.mouse are optional
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor")
cellchat@DB <- CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat)
#future::plan("multiprocess", workers = 24)
# Identify overexpressed genes
cellchat <- identifyOverExpressedGenes(cellchat)
#Identify ligand-receptor pairs
cellchat <- identifyOverExpressedInteractions(cellchat)
#Map ligands and receptors into the PPI network
cellchat <- projectData(cellchat, PPI.human)

## Interaction inference
## 1、calculate the communication probability to infer the communication network of cell interaction
cellchat <- computeCommunProb(cellchat, raw.use = F,population.size = TRUE) 
###Filter out small cell populations
cellchat <- filterCommunication(cellchat, min.cells = 3)

#Extract inferred cell interaction communication network data matrix, which includes all inferred ligand/receptor level cell-cell communication.
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"net_lr_ECM.csv")
levels(cellchat@idents)

##2、Infer intercellular communication at the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) #计算一个通路上的所有配体受体对
##Summarize communication probabilities to calculate aggregated communication networks between cells
cellchat <- aggregateNet(cellchat)

##3、Calculate the aggregation cell interaction communication network
groupSize <- as.numeric(table(cellchat@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE) 
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

##############Save results
pathway.show.all=cellchat@netP$pathways
##Generate bubble plot
levels(cellchat@idents)
dir.create("./bubble")
setwd("./bubble/")
for (i in 1:length(levels(cellchat@idents))) {
  p=netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:length(levels(cellchat@idents))), thresh = 0.05,remove.isolate = F,angle.x = 45)
  ggsave(filename = paste0(gsub("/","_",levels(cellchat@idents)[i]),"_bubble.pdf"),
         p,width=8,height=12,dpi=300)
}

##Individual cells interact with other cells
#a.intensity or probability of interactions
setwd("../")
dir.create("./circle")
setwd("./circle/")
dir.create("./weight")
setwd("./weight/")

mat <- cellchat@net$weight
#par(mfrow = c(1,1), xpd=TRUE)
#par("mar")
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  pdf(file = paste0(gsub("/","_",levels(cellchat@idents)[i]),"_circle_weight.pdf"),width=8,height=12)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,vertex.label.cex = 0.5, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

#b.number of interaction
setwd("../")
dir.create("./count")
setwd("./count/")

mat <- cellchat@net$count
#par(mfrow = c(1,1), xpd=TRUE)
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  pdf(file = paste0(gsub("/","_",levels(cellchat@idents)[i]),"_circle_count.pdf"),width=8,height=12)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,vertex.label.cex = 0.5, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

