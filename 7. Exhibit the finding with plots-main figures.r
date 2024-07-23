##############################Code for Figure 2e######################################################
#load the packages
library(ggplot2)
library(tidyverse)
library(ggrepel)

#load data
df <- read.csv("TF_enrichment_result.csv",header = T)
head(df)

#add a significance label
df$label <- ifelse(df$p_val<0.05,"P-value<0.05","P-value>=0.05")
head(df)

#obtain the significant genes in each cluster
sig <- df[which(df$P<0.05),]
sig<-sig[order(sig$Cancer_ID,-sig$OR),]

#add a new columnto label the significant genes as 2 and the others as 1.
df$size <- case_when(!(df$geneID %in% sig$geneID)~ 1,
df$geneID %in% sig$geneID ~ 2)

#extract the not significant gene dataframe；
dt <- filter(df,size==1)
head(dt)

#draw scatter volcano maps  for not significant genes in each Cluster
p <- ggplot()+
geom_jitter(data = dt,
aes(x = cluster, y = log2FC, color = label),
size = 0.85,
width =0.4)

#Overlay significant gene scatter points of each Cluster, and enlarge these points
p <- ggplot()+
geom_jitter(data = dt,
aes(x = cluster, y = log2FC, color = label),
size = 0.85,
width =0.4)+
geom_jitter(data = sig,
aes(x = cluster, y = log2FC, color = label),
size = 1,
width =0.4)
p

#determine the length of the background column according to the log2OR interval in Figure p：
dfbar<-data.frame(x=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),
                  y=c(3.17,6.38,3.90,4.14,4.65,3.74,2.63,3.81,3.45,3.78,4.00,3.36,6.14,2.57,3.92,2.60,2.28,4.21,1.81,1.83,3.25,4.71,5.06))
dfbar1<-data.frame(x=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),
                   y=c(-0.050 ,-0.050 ,-0.050,-0.05 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050 ,-0.050))

#draw background column
p1 <- ggplot()+
geom_col(data = dfbar,
mapping = aes(x = x,y = y),
fill = "#dcdcdc",alpha = 0.6)+
geom_col(data = dfbar1,
mapping = aes(x = x,y = y),
fill = "#dcdcdc",alpha = 0.6)

#Overlay the scatter volcano map onto the background column
p2 <- ggplot()+
geom_col(data = dfbar,
mapping = aes(x = x,y = y),
fill = "#dcdcdc",alpha = 0.6)+
geom_col(data = dfbar1,
mapping = aes(x = x,y = y),
fill = "#dcdcdc",alpha = 0.6)+
geom_jitter(data = dt,
aes(x = cluster, y = log2FC, color = label),
size = 0.85,
width =0.4)+
geom_jitter(data = sig,
aes(x = cluster, y = log2FC, color = label),
size = 1,
width =0.4)

#Add the cluster color block label for the X axis
dfcol<-data.frame(x=c(0:22),
y=0,
label=c("BLCA","BRCA","CHOL","CRC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SARC","SKCM","STAD","THCA","UCEC","UVM"))
#Assign colors
mycol <- c("#F39696","#CB181D","#A50F15","#d94801","#E2641E","#BE8F4D","#F9C481","#A6C5C9","#528B8B","#385D60","#385B5E","#425A7A","#08519C","#2171B5","#4292C6","#6BAED6","#BCBDDC","#9D9BD1","#807DBA","#6A51A3","#54278F","#3F007D","#2C2C6B")
p3 <- p2 + geom_tile(data = dfcol,
aes(x=x,y=y),
height=0.4,
color = "black",
fill = mycol,
alpha = 0.6,
show.legend = F)

#tag the significant genes in each Cluster
p4 <- p3+
geom_text_repel(
data=sig,
aes(x=cluster,y=log2FC,label=geneID),
force = 1.2,
arrow = arrow(length = unit(0.008, "npc"),
type = "open", ends = "last")
)

#adjust scatter color 
p5 <- p4 +
scale_color_manual(name=NULL,
values = c("red","black"))

#Modify X/Y axis headings and add cluster numbers
p6 <- p5+
labs(x="Cluster",y="average logFC")+
geom_text(data=dfcol,
aes(x=x,y=y,label=label),
size =6,
color ="white")

#Custom theme
p7 <- p6+
theme_minimal()+
theme(
axis.title = element_text(size = 13,
color = "black",
face = "bold"),
axis.line.y = element_line(color = "black",
size = 1.2),
axis.line.x = element_blank(),
axis.text.x = element_blank(),
panel.grid = element_blank(),
legend.position = "top",
legend.direction = "vertical",
legend.justification = c(1,0),
legend.text = element_text(size = 15)
)
p7
ggsave("TF_enrichment_scatter_pot.pdf", p7, height = 5, width = 10, device = "pdf")



##############################Code for Figure 3h######################################################
library(fgsea)
keytypes(org.Hs.eg.db)
#a.hallmark####

#load("pathways.rda")#immune
#Load information of hallmark pathway
load("/home/CYM/data/immunQTL/database/Annotation/gene_correlation/add_165CRC/hallmark/pathways.RData")#hallmark,toType = c("ENTREZID")

#conduct GSEA analysis
k=0.95
fgseaRes_all <- c()
for(a in 1:nrow(cancer_type)){ 
  gene_cor<-read.csv(paste0(cancer_type[a,1],"/","gene_immunqtl_cell_correlation_pcc.csv"),header = T,check.names = F)
  gene<-gene_sig[which(gene_sig$trait==cancer_type[a,1]),]#extract the significante genes
  gene0<-left_join(gene,gene_cor,by=c("Cell.type"="cell_type","symbol"="symbol"))
  gene0<-na.omit(gene0)
  data<-gene0[,c(13,5,17)]
  data<-subset(data,!duplicated(data$symbol) | !duplicated(data$Cell.type))
  cell<-unique(data$Cell.type)

  for(b in 1:length(cell)){
  RS<- data[which(data$Cell.type==cell[b]),]
  if(nrow(RS)<2){next;}
  ids<- bitr(RS$symbol, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db,drop = T)
  ranks<-merge(RS,ids,by.x="symbol",by.y="SYMBOL",sort=F)
  ranks2<-ranks$Rankscore
  names(ranks2)<-ranks$ENTREZID
  
  fgseaRes <- fgsea(pathways =pathways, stats =ranks2, minSize=1, maxSize=5000, nperm=1000) #minSize大于1可能会有报错
  sigValue <- c()
  
  for(j in 1:nrow(fgseaRes)){
    if(nrow(fgseaRes) == 0) {
      next()
    } 
    if(fgseaRes$ES[j]>0){
      sig_ij <- 1 - 2*fgseaRes$pval[j]
    }else{
      sig_ij <- 2*fgseaRes$pval[j] - 1
    }
    sigValue <- c(sigValue,sig_ij)
    
  }
  if(nrow(fgseaRes) == 0) {
    next()
  } 
  trait <- cancer_type[a,1]
  cell_type<-cell[b]
  fgseaRes_i <- cbind(trait,cell_type,fgseaRes,sigValue)
  fgseaRes_all <- rbind(fgseaRes_all,fgseaRes_i)
  }
}
all_pairs <- fgseaRes_all[1:nrow(fgseaRes_all),c(1:9,11)] %>% as.data.frame()
all_pairs_sig<-all_pairs[which(all_pairs$pval<0.05),]

#Start to plot the graph
#select the top20 significant pathway
data<-all_pairs[which(all_pairs$pathway %in% pathway),]
data<-data[order(data$pval,decreasing = F),]
data$ID<-paste(data$trait,data$pathway,sep="_")
data<-data[!duplicated(data$ID),]
data$trait<-factor(data$trait,levels=rev(sort(cancer_type$Cancer)))

p_gsea2 <- ggplot() +
  geom_tile(
    data = data,
    aes(pathway,trait, fill = NES), 
    colour = "white", size = 1
  ) +
  scale_fill_gradient2(low = "#003366",high = "firebrick3",midpoint=0,na.value = "white",guide = "colourbar",
                       name = "NES"
  )+
  geom_point(
    data = data,
    aes(pathway,trait, size = -log10(data$pval)),
    shape = 1
  ) +
  scale_size_area(
    breaks = c(1.5, 2),
    labels = c(1.5, "> 2"),
    max_size=5,
    name = "log10(P)"
  ) + 
  labs(
    x = "",
    y = ""
  ) +
  scale_x_discrete(position = "bottom") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 90),
    legend.position = "right", 
    legend.direction="vertical",
    legend.title = element_text(angle = 90),
    legend.title.align = 0.5,
    legend.box.just = "left"
  ) +
  guides(
    fill = guide_colorbar(title.position = "left", order = 1),
    size = guide_legend(title.position = "left", order = 2)
  )
p_gsea2
ggsave("hallmarker_heatmap_bubble_plot.pdf", p_gsea2, height = 7, width = 8, device = "pdf")

##############################Code for Figure 5k-n######################################################
library(ggalluvial)
#load data
df <- read.table("cell_drug_correlation_pcc.txt",sep = "\t",header = T)
df<-df[,c(2,9,3)]

#assign the colors
mycol <-c("#E5A2A3","#C14D4E","#931F23","#CD6E4A","#D77F51","#BE8F4D","#E8B97E","#A2C0C4","#739F9E","#547F86","#385B5E","#3D5673","#396DA0","#5185B1","#689BBE","#85B1CC","#C5C4DB","#A8A6CC","#918FB9","#7E6DA4","#6F5097","#59517E","#515581")

#Format conversion
UCB_lodes <- to_lodes_form(df[,1:ncol(df)],
                           axes = 1:ncol(df),
                           id = "group")
View(UCB_lodes)

#Start to plot the graph
plot<-ggplot(UCB_lodes,
             aes(x = x, stratum = stratum, alluvium = group,
                 fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) + 
  geom_flow(width = 1/8) + #Set the width of the gap between the line and the square
  geom_stratum(alpha = .9,width = 1/10) + #Set the transparency and width of the block
  geom_text(stat = "stratum", size = 1,color="black",angle=90) + 
  scale_fill_manual(values = mycol) +
  xlab("") + ylab("") +
  theme_bw() + #Remove background color
  theme(panel.grid =element_blank()) + #Remove gridlines
  theme(panel.border = element_blank()) + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + #Remove axis
  ggtitle("")+
  guides(fill = FALSE) 
plot
ggsave(plot,filename ="drug_sankey.pdf",width = 5,height = 6)



##############################Code for Figure 5k-n######################################################
##############rs1360948-survival################
library(survival)
library("readr")
library(stringr)

#load data
geno<-read.table("genotype_data.txt",sep="\t",header = T,check.names = F,row.names = 1)
data <- read.csv("survival_status.csv",header = T,check.names = F,row.names = 1)
data<-t(data)

#extract the part where the sample number overlaps
sample<-colnames(geno)
sample1<-str_split(sample,'_',simplify = T)[,1]
colnames(geno)<-sample1
sample2<-intersect(sample1,colnames(data))
geno<-geno[,sample2]
data<-data[,sample2]

#merge data
x=geno[which(rownames(geno) == "rs1360948"),]
colnames(x)<-colnames(data)
data2=as.data.frame(t(rbind(x,data)))
data2=data2[!is.na(data2[,1]) & !is.na(data2[,2])& !is.na(data2[,3]),]
colnames(data2)[1]=c("snp")
data2=data2[as.numeric(as.character(data2[,3]))>0,]		
gg=c("AA","Aa","aa")
names(gg)=c("0","1","2")

#Convert days to months
data2$futime=as.numeric(as.character(data2[,3]))/30 
data2$fustat=ifelse(as.vector(data2[,2]) == "1", 1,0)
model1 = survdiff(Surv(futime, fustat) ~ snp, data=data2, na.action=na.exclude)
fit1   = survfit(Surv(futime,fustat)~ snp,data=data2,na.action=na.exclude)
p=format(1-pchisq(model1$chisq, df=length(levels(factor(data2$snp)))-1),digits= 3)		

#Start to plot the graph
outfile=c("CRC_survival_plot.pdf")
pdf(outfile,width = 5,height = 5)
plot(fit1,xlab="Months", col=c("#e4a59c","#98cdd8","#7c89a9"),mark.time=TRUE,ylab="Survival Probability",lwd = 1.5,	
     main=paste(Cond,", ",snp,", KM curve" ,"\n","p-value = ",p, sep = ""))

legend("bottomleft", col=c("#e4a59c","#98cdd8","#7c89a9"),  lty=1,bty="n",
       legend= paste(as.vector(gg[attributes(as.factor(data2$snp))$levels]), ", n=", model1$n, sep = ""))

dev.off()


##############################Code for Figure 7c######################################################
#conduct survival analysis
covariance<-read.csv("TCGA_CRC_survival_data.csv",header = T,check.names = F,row.names = 1)
library(survival)
my.surv <- Surv(covariance$days_to_last_follow_up/30,covariance$`Survival status`)
#list the covariates
age<-covariance$Age 
gender<-covariance$Gender
stage<-covariance$Stage 

survival_dat<-data.frame(values1=covariance$group, gender=gender,stage=stage,age=age, stringsAsFactors = F) #
m=coxph(my.surv~age+gender+stage+values1,data=survival_dat) #
beta<-coef(m)
se<-sqrt(diag(vcov(m)))
HR<-exp(beta)
HRse <- HR * se
survival_result <- signif(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                          HR = HR, HRse = HRse,
                          HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                          HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                          HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)#0.298

#start to plot the graph
fit=survfit(my.surv~values1,data=survival_dat) 
ggsurv <- ggsurvplot(
  fit,                     
  #risk.table = TRUE,       
  pval = survival_result$p[4],             
  conf.int = TRUE,         
  palette = c("#3A669A","#B44D4D"),
  #xlim = c(0,150),         
  xlab = "Time in months",   
  break.time.by = 50,     
  ggtheme = theme_light(), 
  risk.table.y.text.col = T,
  risk.table.height = 0.25, 
  risk.table.y.text = FALSE,
  #ncensor.plot = TRUE,      
  #ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  
  legend.labs = c("Low", "High")    
)
ggsurv
