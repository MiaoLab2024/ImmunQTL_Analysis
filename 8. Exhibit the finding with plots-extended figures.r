##############################Code for extended Figure 6a######################################################
#load data
data_GO$cancer_type<-factor(data_GO$cancer_type,levels=rev(sort(cancer_type$Cancer)))  
data_GO$Description<-factor(data_GO$Description,levels=GO_sig$Description)  

#Start to plot the graph
g3_2 =ggplot(data_GO, aes(x =Description ,  y = cancer_type)) + 
  geom_point(aes(size = Count,fill=-log10(data_GO$pvalue)), shape = 21,stroke=0) + 
  scale_size_continuous(range=c(3,10))+ 
  scale_fill_gradient2(low = "white",high = "#e7535f",limits=c(0,5))+ #4b699d blue; #6dc361 GREEN;#e7535f red
  theme_minimal()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray",size = rel(0.5)),
        panel.border = element_rect(colour="black",fill=NA),
        axis.title.x=element_blank(),#remove title
        axis.ticks.x=element_blank(),#remove x-axis
        axis.title.y=element_blank(),#remove y-axis
        axis.text.x = element_text(angle = 90, hjust = 1)) # 调整x轴文字 )
g3_2 
ggsave(g3_2 ,filename = "ALL_GO_BP_bubble_plot.pdf",  width = 8, height = 8)



##############################Code for extended Figure 6c######################################################
all_cell_survival<-read.csv("all_cell_survival_cox_for_plot.csv",header = T,check.names = F)
factor<-rev(cancer_type[,1])
all_cell_survival$cancer_type = factor(all_cell_survival$cancer_type, levels=factor)

#Start to plot the graph
library(RColorBrewer)
cols<-c("red"="#EA5D66","blue"="#546FA2")
ggplot(cell_survival) +
  geom_point(aes(x = cell_survival$ID,y=cell_survival$cancer_type,size=cell_survival$prop,fill = cell_survival$color), pch = 21,color = "white",alpha=0.9) +
  scale_fill_manual(values = cols,name="HR") +
  scale_size(range = c(2,8),name="Prop")+
  scale_x_discrete(expand = expand_scale(mult = c(0.02,0.02)),name = "Cell types")+
  theme_bw() +
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.text.x = element_text(angle = 45,hjust=1),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 8, face = "bold"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    legend.background = element_blank(),
    panel.grid.major = element_line(size = .5),
  )

ggsave("cell_survival_bubbel_plot.pdf", width = 8, height = 7)    


