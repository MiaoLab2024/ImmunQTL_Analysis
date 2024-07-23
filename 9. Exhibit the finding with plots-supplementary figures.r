##############################Code for Fig.S3c######################################################
#load the packages
library(ggplot2)
library(ggrepel)

#load data
result<-read.csv("rbp_enrichment_result.csv",header = T,check.names = F,row.names = 1)
level<-cancer_type$Cancer
result$Cancer=factor(result$Cancer,levels=rev(level))
result$log10P<-(-log10(result$P))
result$log2OR<-log2(result$OR)

#Start to plot the graph
p <- ggplot() +
  geom_tile(
    data = result,
    aes(Experiment.target, Cancer, fill = log2OR), 
    colour = "white", size = 1
  ) +
  scale_fill_gradient2(low = "#003366",high = "firebrick3",midpoint=0,na.value = "white",guide = "colourbar",
                       name = "log2(OR)"
  )+
  geom_point(
    data = result,
    aes(Experiment.target, Cancer, size = log10P),
    shape = 1
  ) +
  scale_size_area(
    breaks = c(1, 2, 3),
    labels = c(1, 1.5, ">2"),
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

p
ggsave("RBP_heatmap_bubble_plot.pdf", p, height = 5, width = 12, device = "pdf")
