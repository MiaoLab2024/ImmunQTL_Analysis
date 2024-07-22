library(nichenetr)
##Load prior data for nichenet analysis
lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")

dir.create("./ligand-target")
dir.create("./ligand-receptor")

Idents(sce.ESCA) <- "Immune_All_Low"
cell<-c("Cycling B cells","DC1","Endothelial cells","Epithelial cells","Fibroblasts","gamma-delta T cells","Memory B cells","NK cells","Tcm/Naive cytotoxic T cells","Tcm/Naive helper T cells","Type 17 helper T cells","Macrophages","Plasma cells")

for(i in 1:length(cell)){
  nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = sce.ESCA, 
                                                 top_n_ligands = 20,
                                                 receiver = cell[i], 
                                                 sender = "all",
                                                 condition_colname = "orig.ident", 
                                                 condition_oi = "Tumor", 
                                                 condition_reference = "Normal", 
                                                 ligand_target_matrix = ligand_target_matrix, 
                                                 lr_network = lr_network, 
                                                 weighted_networks = weighted_networks)
  
  ## View the ligand activity analysis results
  # Mainly refer to pearson index, bona_fide_ligand=True represents the ligand-receptor reported in the literature.
  # bona_fide_ligand=False means that PPI predicts ligand-receptors that have not been experimentally proven.
  x <- nichenet_output$ligand_activities
  ligand_aupr_matrix = as.matrix(x$aupr_corrected) 
  rownames(ligand_aupr_matrix)<-x$test_ligand
  order_ligands = rownames(nichenet_output$ligand_target_matrix)
  order_ligands <- gsub("\\.", "-", order_ligands)
  #order_ligands<-c("TNFSF13","HLA-DQA2") i=9
  vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
  p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)")
  #p_ligand_aupr
  x$receiver<-cell[i]
  write.csv(x, paste0(gsub("/","_",cell[i]),"_ligand_activities.csv"), row.names = F)
  p_ligand_expression<-nichenet_output$ligand_expression_dotplot
  p_ligand_expression<-p_ligand_expression+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0(gsub("/","_",cell[i]),"_Dotplot_ligand-expression.pdf"), p_ligand_expression, width = 12, height = 6)
  
  
  ## View ligand-regulated target genes
  p_ligand_target = nichenet_output$ligand_target_heatmap+ 
  scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + 
  ggsave(paste0(gsub("/","_",cell[i]),"_Heatmap_ligand-target.pdf", p, width = 12, height = 6)
  # View the top target genes regulated by ligand and their scores
  x <- nichenet_output$ligand_target_matrix
  write.csv(x, paste0("ligand-target/",gsub("/","_",cell[i]),"_ligand_target_matrix.csv"), row.names = F)
  
  ##Combine all figures
  figures_without_legend=cowplot::plot_grid(
    p_ligand_aupr+theme(legend.position="none",axis.ticks=element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
    p_ligand_target + theme(legend.position = "none",axis.ticks = element_blank ()) + ylab(""),
    align ="hv",nrow = 1
  )
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
    nrow = 1,
    align = "h")
  
  combine_plot<-cowplot::plot_grid(figures_without_legend, 
                                   legends, 
                                   rel_heights = c(10,2), nrow = 2, align = "hv")
  ggsave(paste0("ligand-target/",gsub("/","_",cell[i]),"_Heatmap_ligand-target.pdf"), combine_plot, width = 12, height = 6)
  
  
  ## View the receptors
  # View ligand-receptor interactions
  p_ligand_receptor = nichenet_output$ligand_receptor_heatmap
  ggsave(paste0("ligand-receptor/",gsub("/","_",cell[i]),"_Heatmap_ligand-receptor.pdf"), p_ligand_receptor, width = 12, height = 6)
  x <- nichenet_output$ligand_receptor_matrix
  write.csv(x, paste0("ligand-receptor/",gsub("/","_",cell[i]),"_ligand_receptor_matrix.csv"), row.names = F)
  
  # View the ligand-receptor reported in the literature
  # Show ‘bona fide’ ligand-receptor links that are described in the literature and not predicted based on PPI
  p = nichenet_output$ligand_receptor_heatmap_bonafide
  #ggsave("DC2_Heatmap_ligand-receptor_bonafide.pdf", p, width = 8, height = 4)
  x <- nichenet_output$ligand_receptor_matrix_bonafide
  #x <- nichenet_output$ligand_receptor_df_bonafide
  write.csv(x, paste0("ligand-receptor/",gsub("/","_",cell[i]),"_ligand_receptor_bonafide.csv"), row.names = F)
}


