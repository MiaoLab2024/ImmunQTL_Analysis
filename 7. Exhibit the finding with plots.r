################################Code for figure 2e##########################################################
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

#obtain the top10 significant genes in each cluster
top10sig0 <- filter(df,cluster=="0") %>% distinct(geneID,.keep_all = T) %>% top_n(10,abs(log2FC))
head(top10sig0)

#Combine the top10 gene dataframes from all clusters
top10sig <- rbind(top10sig0,top10sig1,top10sig2,top10sig3,top10sig4,top10sig5,top10sig6,top10sig7,top10sig8)

#add a new columnto label the Top10 genes as 2 and the others as 1.
df$size <- case_when(!(df$geneID %in% top10sig$geneID)~ 1,
df$geneID %in% top10sig$geneID ~ 2)

#extract the non-top10 gene dataframe；
dt <- filter(df,size==1)
head(dt)

#draw scatter volcano maps  for non-top10 genes in each Cluster
p <- ggplot()+
geom_jitter(data = dt,
aes(x = cluster, y = log2FC, color = label),
size = 0.85,
width =0.4)

#Overlay Top10 gene scatter points of each Cluster, and enlarge these points
p <- ggplot()+
geom_jitter(data = dt,
aes(x = cluster, y = log2FC, color = label),
size = 0.85,
width =0.4)+
geom_jitter(data = top10sig,
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
geom_jitter(data = top10sig,
aes(x = cluster, y = log2FC, color = label),
size = 1,
width =0.4)

#Add the cluster color block label for the X axis
dfcol<-data.frame(x=c(0:22),
y=0,
label=c("BLCA","BRCA","CHOL","CRC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SARC","SKCM","STAD","THCA","UCEC","UVM"))
mycol <- c("#F39696","#CB181D","#A50F15","#d94801","#E2641E","#BE8F4D","#F9C481","#A6C5C9","#528B8B","#385D60","#385B5E","#425A7A","#08519C","#2171B5","#4292C6","#6BAED6","#BCBDDC","#9D9BD1","#807DBA","#6A51A3","#54278F","#3F007D","#2C2C6B")
p3 <- p2 + geom_tile(data = dfcol,
aes(x=x,y=y),
height=0.4,
color = "black",
fill = mycol,
alpha = 0.6,
show.legend = F)

#tag the Top10 genes in each Cluster
p4 <- p3+
geom_text_repel(
data=top10sig,
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
