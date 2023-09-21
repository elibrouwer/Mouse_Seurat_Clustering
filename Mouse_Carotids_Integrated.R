
# Mouse sample corrected for tissue  --------------------------------------

Human_Mouse_Obj.combined.list <- SplitObject(Human_Mouse_Obj.combined, split.by = "dataset")

Mouse_tissue_cor <- Human_Mouse_Obj.combined.list[[1]]

Mouse_tissue_cor@meta.data %<>%
  dplyr::mutate(SeuratIntegrated = JakeLusis$PlaqueCormat,
                Batch = Strain_Batch$Batch,
                Strain = JakeLusis$Strain)  

Mouse_tissue_cor %<>% #logtransformation natural log
  ScaleData( verbose = FALSE, features = rownames(Mouse_Obj.combined)) %>%
  RunPCA( verbose = FALSE)  %>%
  FindNeighbors( dims = 1:12) %>%
  FindClusters( resolution = 0.448) %>%
  RunUMAP( dims = 1:12) %>%
  RunTSNE( dims = 1:12)

MouseUmap_tis <- DimPlot(Mouse_tissue_cor, reduction = "umap", group.by = "seurat_clusters", pt.size = 1.5) + labs(title = "Mouse tissue corrected Dataset clustered with Seurat")
MouseUmap_tis

MouseUmap_Plaque_tis <- DimPlot(Mouse_tissue_cor, reduction = "umap", group.by = "SeuratIntegrated", pt.size = 1.5) + labs(title = "Mouse tissue corrected Dataset clustered with Seurat")
MouseUmap_Plaque_tis

MouseUmap_tis + MouseUmap_Plaque_tis

Mouse_tissue_cor@meta.data

FindAllMarkers(Mouse_tissue_cor, only.pos= TRUE, min.pct = 0.2, logfc.threshold = 0.1) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)  -> top10_tis

DiffExpHeat_tis <- DoHeatmap(Mouse_tissue_cor, features = top10_tis$gene)  + scale_fill_gradientn( colors = colorpanel(100,"grey100","grey90","darkgreen")) + labs(title = "tissue corrected Differential Gene Expression Mouse Dataset, MouseClusters")
DiffExpHeat_tis

DiffExpHeat_tis_HumanPlaq <- DoHeatmap(Mouse_tissue_cor, features = top10_tis$gene, group.by = "SeuratIntegrated")  + scale_fill_gradientn( colors = colorpanel(100,"grey100","grey90","darkgreen")) + labs(title = "Differential Gene Expression Mouse Dataset, MouseClusters")
DiffExpHeat_tis_HumanPlaq

DiffExpHeatFig1_tis <- DoHeatmap(Mouse_tissue_cor, features = figure1genelist$GenesSymbol)  + scale_fill_distiller(palette = "Greens", direction = 1) + labs(title = "Differential Gene Expression Mouse Dataset with important Human genes")
DiffExpHeatFig1_tis

DiffExpHeatFig1_tis_HumanPlaq <- DoHeatmap(Mouse_tissue_cor, features = figure1genelist$GenesSymbol, group.by = "SeuratIntegrated")  + scale_fill_gradientn( colors = colorpanel(100,"grey100","grey90","darkgreen")) + labs(title = "Differential Gene Expression Mouse Dataset, MouseClusters")
DiffExpHeatFig1_tis_HumanPlaq
