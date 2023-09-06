#Seurat integration method
if (!require("pacman")) install.packages("pacman")
update.packages("pacman")

pacman::p_load(purrr,gtsummary,scclusteval,gtable,gplots,grid,tidyHeatmap,org.Mm.eg.db,DESeq2,ReactomeContentService4R)
library(purrr)
library(gtsummary)
library(scclusteval)
library(gtable)
library(gplots)
library(grid)
library(tidyHeatmap)
library(org.Mm.eg.db)
library(DESeq2)
library(ReactomeContentService4R)
library(ReactomePA)
library(limma)
library(ggpubr)
library(BiocManager)
library(biomaRt)
library(plyr)
library(dplyr)
library(ROCR)
library(RColorBrewer)
library(ggplot2)
library(magrittr)
library(conflicted)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(DESeq2)
library(tools)
library(gridExtra)
library(stringr)
library(cowplot)
library(harmony)
library(Matrix)
library(ReactomePA)
library(org.Hs.eg.db)
library(biomaRt)
library(Seurat)
library(readr)
library(GenomeInfoDbData)

update.packages()
###Tables that form checkpoints for this file###
# JakelusisTable_With_SeuratClustering ------------------------------------

JakeLusis <- read_csv("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569_JakeLusis_Table_SeurataClustering.csv")
# SeuratIntegration -------------------------------------------------------


MousePath <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569.txt"
HumanPath <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Human/Human1/bulk_RNAseq_raw_counts.txt.minRib.txt.PC.txt"
#Human_Mouse_Obj <- CreateSeuratObject(counts = Mouse_HumanTable(HumanPath, MousePath, Ens_Mouse, Ens_Human)  , project = paste0(Humanversion, "_", DataName))

# Human_Mice Batch Corrected table ----------------------------------------------
Integrated_Mouse <- read_csv("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569_Corrected_Batch.csv") %>%
  column_to_rownames(var = "gene_id")

Human_GeneSymb <- read_delim(file = paste0(file_path_sans_ext(HumanPath), "_", basename(dirname(HumanPath)),"Gene_Symbol" ,".txt")) 
topo.colors(n = 3)
Integrated_Mouse %>%
  rownames_to_column(var = "gene_id") %>%
  left_join(Human_GeneSymb, by = c("gene_id")) %>%
  distinct(gene_id , .keep_all = TRUE) %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'gene_id') %>%
  normalizeQuantiles()  -> Human_Mice

# Human_Mouse integration -------------------------------------------------------------
#Human_Mice <- read_delim(file = paste0(file_path_sans_ext(HumanPath), "_", basename(dirname(HumanPath)),"Human_Mice" ,".txt")) %>%
  column_to_rownames(var = 'gene_id')

Human_Mouse_Obj <- CreateSeuratObject(counts = Human_Mice, project = paste0(Humanversion, "_", DataName))

cluster0 <- SeuratClusters[which(SeuratClusters$cluster == 0), ] #list of human clusters with human ids
cluster1 <- SeuratClusters[which(SeuratClusters$cluster == 1), ]
cluster2 <- SeuratClusters[which(SeuratClusters$cluster == 2), ]
cluster3 <- SeuratClusters[which(SeuratClusters$cluster == 3), ]
cluster4 <- SeuratClusters[which(SeuratClusters$cluster == 4), ]

Human_Mouse_Obj %<>%
  ScaleData( verbose = FALSE, features = rownames(Human_Mouse_Obj)) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = nrow(Human_Mouse_Obj)) %>%
  RunPCA( npcs = 30, verbose = FALSE) %>% 
  RunUMAP( dims = 1:20)


Human_Mouse_Obj@meta.data %<>%
  rownames_to_column("Names") %>%
  dplyr::mutate(dataset = ifelse(grepl("ae",Names), paste0("Human_Caro"), 
                                 ifelse(grepl("GSM",Names), paste0("Mouse_APOE"), "-"))) %>% #add dataset to metadata
  dplyr::mutate(Plaque_Cat = ifelse(Names %in% cluster0$study_number, 0,
                             ifelse(Names %in% cluster1$study_number,1,
                                    ifelse(Names %in% cluster2$study_number, 2,
                                           ifelse(Names %in% cluster3$study_number, 3,
                                                  ifelse(Names %in% cluster4$study_number, 4, NA)))))) %>% #add cluster to metadata
  column_to_rownames("Names")


Human_Mouse_Obj.list <- SplitObject(Human_Mouse_Obj, split.by = "dataset")

Human_Mouse_Obj.list <- lapply(X = Human_Mouse_Obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000000)
})

features <- SelectIntegrationFeatures(object.list = Human_Mouse_Obj.list)

Human_Mouse_Obj.anchors <- FindIntegrationAnchors(object.list = Human_Mouse_Obj.list, anchor.features = features)
Human_Mouse_Obj.combined <- IntegrateData(anchorset = Human_Mouse_Obj.anchors)

DefaultAssay(Human_Mouse_Obj.combined) <- "integrated"

Human_Mouse_Obj.combined %<>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>% 
  RunUMAP( dims = 1:12) %>%
  RunTSNE(dims = 1:12)


# readRDS_Human_Obj -------------------------------------------------------
##run object from here 
#last save 8-5
#saveRDS(Human_Mouse_Obj.combined, file = "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/SeuratObjects/HumanJake_Integrated_BatchCorrect.rds")
Human_Mouse_Obj.combined <- readRDS("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/SeuratObjects/HumanJake_Integrated_BatchCorrect.rds")
# eigenvalues_and_clust ---------------------------------------------------

mat <- Seurat::GetAssayData(Human_Mouse_Obj.combined, assay = "integrated", slot = "scale.data")
pca <- Human_Mouse_Obj.combined[["pca"]]

total_variance <- sum(matrixStats::rowVars(mat)) 

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
sum(varExplained[1:20])

# Dimplots ----------------------------------------------------------------
pca_After_Integration <- DimPlot(object = Human_Mouse_Obj.combined, reduction = "pca",pt.size = 1.5,shape.by = c('dataset') , group.by = c( "Plaque_Cat") )

Human_Mouse_Obj@meta.data %<>% 
  dplyr::mutate(Plaque_Cat_Mouse  = replace_na(as.character(Plaque_Cat), "Mouse_APOE")) #changing NA's into Mouse_Apoe for the plots
Human_Mouse_Obj.combined@meta.data %<>% 
  dplyr::mutate(Plaque_Cat_Mouse  = replace_na(as.character(Plaque_Cat), "Mouse_APOE")) #changing NA's into Mouse_Apoe for the plots

APOE_Before_Integration <-  DimPlot(Human_Mouse_Obj, reduction = "umap",shape.by = c('dataset') , group.by = c( "Plaque_Cat_Mouse") , pt.size = 1.5) + labs(title = "UMAP before Seurat integration")
APOE_After_Integration <- DimPlot(Human_Mouse_Obj.combined, reduction = "umap",shape.by = c('dataset') , group.by = c( "Plaque_Cat_Mouse") , pt.size = 1.5) + labs(title = "UMAP after Seurat integration")

APOE_Before_Integration + APOE_After_Integration
APOE_After_Integration_Tsne <- DimPlot(Human_Mouse_Obj.combined, reduction = "tsne",shape.by = c('dataset') , group.by = c( "Plaque_Cat_Mouse") , pt.size = 1.5) + labs(title = "tnse after Seurat integration")

Mouse_control <- Human_Mouse_Obj@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column(var = "ID") %>%
  dplyr::filter(str_detect(ID, "^GSM")) %>%
  left_join(JakeLusis, by = c("ID" = "categories"))

Human_Mouse_Obj.combined@meta.data %<>% #adding the new col "corrected" Plaque Cat, this is the plaque cat from the pearsons correlation
  rownames_to_column(var = 'categories') %>%
  left_join(JakeLusis , by = 'categories') %>%
  dplyr::select(-c(cluster,comparison,Strain,Batch, value, SeuratCluster, SeuratCluster.x, SeuratCluster.y)) %>%
  dplyr::mutate(Corrected_Plaque_Cat = coalesce(Plaque_Cat, PlaqueCormat)) %>%
  column_to_rownames(var = 'categories') 

APOE_After_Correlation <- DimPlot(Human_Mouse_Obj.combined, reduction = "umap",shape.by = c('dataset') , group.by = c( "Corrected_Plaque_Cat") , pt.size = 1.5) + labs(title = "UMAP after Seurat integration")
# Correlate Mice to human -------------------------------------------------------------

APOE_cormat <- as.data.frame(stats::cor(as.matrix(Seurat::GetAssayData(Human_Mouse_Obj.combined, assay = "integrated", slot = "data")))) %>% #pearsons correlation of Human and APOE Mouse
  tibble::rownames_to_column(var = "Names") 

APOE_Mice <- tibble(.rows = length(APOE_cormat)-1)

APOE_Mice %<>%
  dplyr::mutate(Names = APOE_cormat$Names) %>%
  left_join(APOE_cormat, by = c("Names")) %>%
  dplyr::mutate(
    Plaque_Cat = ifelse(Names %in% cluster0$study_number, 0, 
                        ifelse(Names %in% cluster1$study_number,1, 
                               ifelse(Names %in% cluster2$study_number, 2, 
                                      ifelse(Names %in% cluster3$study_number, 3, 
                                             ifelse(Names %in% cluster4$study_number, 4, NA)))))) %>% #Adding the plaque cluster from human clustering
  na.omit() 

APOE_Mice %>%
  tidyr::pivot_longer(!c(Plaque_Cat, Names), names_to = "categories", values_to = "Correlation") %>%
  group_by(categories) %>%
  tidyr:: nest() %>%
  dplyr::mutate(cor_plaque0 = map(data, ~dplyr::filter(.x, Plaque_Cat == 0)),
                cor_plaque1 = map(data, ~dplyr::filter(.x, Plaque_Cat == 1)),
                cor_plaque2 = map(data, ~dplyr::filter(.x, Plaque_Cat == 2)),
                cor_plaque3 = map(data, ~dplyr::filter(.x, Plaque_Cat == 3)),
                cor_plaque4 = map(data, ~dplyr::filter(.x, Plaque_Cat == 4)),) %>% 
  dplyr::mutate("0" = map_dbl(cor_plaque0, ~median(.x$Correlation)),
                "1" = map_dbl(cor_plaque1, ~median(.x$Correlation)),
                "2" = map_dbl(cor_plaque2, ~median(.x$Correlation)), 
                "3" = map_dbl(cor_plaque3, ~median(.x$Correlation)),
                "4" = map_dbl(cor_plaque4, ~median(.x$Correlation))) %>% #adding the correlation to every different plaque type
  dplyr::select(-c(cor_plaque0, cor_plaque1, cor_plaque2, cor_plaque3, cor_plaque4)) %>%
  pivot_longer(-c(categories, data), names_to = "PlaqueCormat") %>%
  group_by(categories) %>% 
  slice_max(value) %>% 
  ungroup() %>%
  left_join(SeuratClusters, by = c( "categories" = "study_number" )) %>%
  dplyr::mutate(comparison = ifelse(PlaqueCormat > cluster, "X",
                                    ifelse(PlaqueCormat < cluster, "X", ""))) %>% #step to control if the boxplots were correct
  dplyr::select(-c(data)) -> Mice_CormatPlaque

APOE_plots <- APOE_Mice %>%
  tidyr::pivot_longer(!c(Plaque_Cat, Names), names_to = "categories", values_to = "Correlation") %>%
  group_by(categories) %>%
  tidyr:: nest() %>%
  dplyr::mutate(boxplot = purrr::map2(data, categories, ~ ggplot(data = .x , aes(x = factor(Plaque_Cat), y = Correlation, fill = factor(Plaque_Cat)))
                                      + geom_boxplot()+ xlab(.y)
                                      + ylab("Correlation")  
                                      + scale_y_continuous(limits = c(0.0, 1))
                                      + scale_fill_brewer(palette = "Paired") 
                                      + theme(legend.position = "none")  
                                      + stat_compare_means( method = "anova", label.y = .95, label.y.npc = "top", label.x.npc = "middle") 
                                      + stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", label.y = .90) )) %>%
  left_join(Mice_CormatPlaque, by = "categories") %>%
  dplyr::mutate(annotateBox = purrr::map2(comparison, boxplot, ~.y+ annotate(geom = "label", label = paste0("Highest median not equal to SeuratCluster : ", comparison), x = 2, y = .3, fill="white"))) %>%
  dplyr::mutate(annotateBox = purrr::map2(cluster, annotateBox, ~.y+ annotate(geom = "label", label = paste0("SeuratCluster Identified ", cluster), x = 2, y = .4, fill="white")))

JakeLusis <- APOE_plots[0:366,] #selecting only the mouse samples 

marranged_Plots <- marrangeGrob(APOE_plots$boxplot[0:367],  nrow = 6, ncol = 8 ) #saving correlation boxplots

ggsave(paste0(DataName, "Mouse SeuratIntegration.pdf"), marranged_Plots, width = 100, height = 80, units = "cm") 

# StrainMice_PlaqueType ---------------------------------------------------
tableJake <- table(JakeLusis$PlaqueCormat)
tableJake

JakeLusis %>%
  dplyr::select(categories, Strain) %>% 
  dplyr::group_by(Strain) %>%
  nest() %>%
  dplyr::mutate(n_total = map(data, nrow)) %>%
  dplyr::select(-c(data)) -> n_strain

mouse_plaque4 <- JakeLusis %>%
  dplyr::filter(PlaqueCormat == "4") %>%
  dplyr::select(categories, Strain) %>% 
  dplyr::group_by(Strain) %>%
  nest() %>%
  dplyr::mutate(n_plaque4 = map(data, nrow)) %>%
  left_join(n_strain, by = c("Strain")) %>%
  dplyr::select(-c(data)) %>%
  as.data.frame()

mouse_plaque4 <- lapply(mouse_plaque4, unlist)
str(mouse_plaque4)
mouse_plaque1 <- JakeLusis %>%
  dplyr::filter(PlaqueCormat == "1") %>%
  dplyr::select(categories, Strain) %>%
  dplyr::group_by(Strain) %>%
  nest() %>%
  dplyr::mutate(n_plaque1 = map(data, nrow)) %>%
  left_join(n_strain, by = c("Strain")) %>%
  dplyr::select(-c(data))
mouse_plaque1 <- lapply(mouse_plaque1, unlist)

#write.table(mouse_plaque4, "MouseStrains_HumanPlaque4.csv", row.names = F)
#write.table(mouse_plaque1, "MouseStrains_HumanPlaque1.csv", row.names = F)


# Analyis of batch corrected mice and human object ------------------------
##run object from here after batch correction
Mouse_Obj.combined <- readRDS("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/SeuratObjects/Jake_BatchCorrect.rds")


JakeLusis %<>%
  left_join(Strain_Batch, by = c("categories" = "Name")) 

Mouse_Obj.combined@meta.data %<>%
  dplyr::mutate(SeuratIntegrated = JakeLusis$PlaqueCormat,
                Batch = JakeLusis$Batch,
                Strain = JakeLusis$Strain)  #adding batch and strain to Mouse_Obj.Combined


Mouse_Obj.combined %<>% #logtransformation natural log
  ScaleData( verbose = FALSE, features = rownames(Mouse_Obj.combined)) %>%
  RunPCA( verbose = FALSE)  %>%
  FindNeighbors( dims = 1:12) %>%
  FindClusters( resolution = 0.448) %>%
  RunUMAP( dims = 1:12) %>%
  RunTSNE( dims = 1:12)

Idents(Mouse_Obj.combined)
RenameIdents(Mouse_Obj.combined, '0' = 'A', '1' = 'B', '2' = 'C', '3' = 'D' ) -> Mouse_Obj.combined
Mouse_Obj.combined@meta.data$seurat_clusters <- Idents(Mouse_Obj.combined) #renaming clusters
# Control_Clustering_FindClusters_DifferentRes -----------------------------------------------
res_list = c(0.1,0.2,0.3,0.4,0.5, 0.448)

for (i in 1:length(res_list)){
  Mouse_Obj.combined %<>%
    FindClusters( resolution = res_list[i])
}

# Calculate_Eigenvalues ---------------------------------------------------

mat <- Seurat::GetAssayData(Mouse_Obj.combined, assay = "integrated", slot = "scale.data")
pca <- Mouse_Obj.combined[["pca"]]

total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  
varExplained = eigValues / total_variance
sum(varExplained[1:12])

Elbowplot_Mice <- ElbowPlot(Mouse_Obj.combined) + labs(title = "ElbowPlot for GSE66569")
Elbowplot_Mice
# Clustering_different_res ------------------------------------------------
reductions_df <- data.frame(row.names = c("pca", "umap", "tsne")) %>%
  rownames_to_column(var = "reduction") %>%
  dplyr::mutate( p3 = map(reduction, ~DimPlot(Mouse_Obj.combined, reduction = paste(.x), label = TRUE, group.by = "integrated_snn_res.0.448", pt.size = 1.5) + ggtitle("res 0.4") + theme(text = element_text(size = 10, face = "bold"))),
                 p1 = map(reduction, ~DimPlot(Mouse_Obj.combined, reduction = paste(.x), label = TRUE, group.by = "integrated_snn_res.0.2", pt.size = 1.5) + ggtitle("res 0.2")),
                 p2 = map(reduction, ~DimPlot(Mouse_Obj.combined, reduction = paste(.x), label = TRUE, group.by = "integrated_snn_res.0.3", pt.size = 1.5) + ggtitle("res 0.3")),
                 p4 = map(reduction, ~DimPlot(Mouse_Obj.combined, reduction = paste(.x), group.by = "integrated_snn_res.0.5", label = TRUE, pt.size = 1.5) +  ggtitle("res 0.5"))) %>%
  pivot_longer(cols = !reduction, values_to = "data") %>%
  group_by(reduction) %>%
  nest() %>%
  dplyr::mutate(data = map(data, ~ .x %>% pivot_wider(names_from = name, values_from = data))) %>%
  dplyr::mutate(res_plots = map2(.x = data, .y = reduction, ~arrangeGrob(grobs = c(.x$p1 , .x$p2 , .x$p3 , .x$p4), nrow = 2, ncol = 2, top = textGrob( paste("Seurat Clustering different resolutions", .y), gp = gpar(fontsize = 20, fontface = 'bold'))))) #clustering at different resolutions

# Silhouette_scores -------------------------------------------------------
Idents(Mouse_Obj.combined)<- Mouse_Obj.combined@meta.data$integrated_snn_res.0.2
silhouette_scores<- CalculateSilhouette(Mouse_Obj.combined, dims = 1:12)
sil_p2<- scclusteval::SilhouetteRainCloudPlot(silhouette_scores)  +  ggtitle("res 0.2")

Idents(Mouse_Obj.combined)<- Mouse_Obj.combined@meta.data$integrated_snn_res.0.448
silhouette_scores<- CalculateSilhouette(Mouse_Obj.combined, dims = 1:12)
sil_p3<- scclusteval::SilhouetteRainCloudPlot(silhouette_scores) +  ggtitle("res 0.448")

Idents(Mouse_Obj.combined)<- Mouse_Obj.combined@meta.data$integrated_snn_res.0.5
silhouette_scores<- CalculateSilhouette(Mouse_Obj.combined, dims = 1:12)
sil_p4<- scclusteval::SilhouetteRainCloudPlot(silhouette_scores) +  ggtitle("res 0.5")

sil_p2 + sil_p3 + sil_p4
# Clustering_plots --------------------------------------------------------

MousePca <- DimPlot(Mouse_Obj.combined, reduction = "pca", group.by = "seurat_clusters")
MousePca

MouseUmap <- DimPlot(Mouse_Obj.combined, reduction = "umap", group.by = "seurat_clusters", pt.size = 1.5) + labs(title = "Mouse Dataset clustered with Seurat")
MouseUmap

MouseUmapSeuratClus <- DimPlot(Mouse_Obj.combined, reduction = "umap", group.by = "SeuratIntegrated", pt.size = 1.5) + labs(title = "Closest correlating Human plaque type")
MouseUmapSeuratClus

MousetsneSeuratClus <- DimPlot(Mouse_Obj.combined, reduction = "tsne", group.by = "seurat_clusters")
MousetsneSeuratClus

MousePca + MouseUmap + MousetsneSeuratClus

MousetsneHumanPlaq <- DimPlot(Mouse_Obj.combined, reduction = "tsne", group.by = "SeuratIntegrated")
MousetsneHumanPlaq



Mouse_Umap_Batch <- DimPlot(Mouse_Obj.combined, reduction = "umap", group.by = "Batch", pt.size = 1.5) + NoLegend()  + labs(title = "Batch After")

Mouse_PCA_Batch <- DimPlot(Mouse_Obj.combined, reduction = "pca", group.by = "Batch", pt.size = 1.5)  + NoLegend()  + labs(title = "Batch After")

Mouse_Umap_Batch + Umap_Before_Batch 

Mouse_PCA_Batch +PCA_Before_Batch 


MouseUmap + MouseUmapSeuratClus

MouseUmapStrain <- DimPlot(Mouse_Obj.combined, reduction = "umap", group.by = "Strain", label = FALSE) 

MouseUmap + MouseUmapStrain


Batch_table <- table(Mouse_Obj.combined@meta.data$seurat_clusters, Mouse_Obj.combined@meta.data$Batch)
Batch_table
# Differential_Expression -------------------------------------------------

FindAllMarkers(Mouse_Obj.combined, assay = "integrated", only.pos= TRUE, min.pct = 0.2, logfc.threshold = 0.1) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)  -> top10 


DiffExpHeat <- DoHeatmap(Mouse_Obj.combined, features = top10$gene)  + scale_fill_gradientn( colors = colorpanel(100,"grey100","grey90","darkgreen")) + plot_annotation(title = "Differential Gene Expression Mouse Dataset", theme = theme(plot.title = element_text(face="bold", size = 18)))
DiffExpHeat #differential expression per MouseCluster

DiffExpHeat_HumanPlaques <- DoHeatmap(Mouse_Obj.combined, features = top10$gene, group.by = "SeuratIntegrated")  + scale_fill_gradientn( colors = colorpanel(100,"grey100","grey90","darkgreen")) + labs(title = "Differential Gene Expression Mouse Dataset, Human Plaque types!")
DiffExpHeat_HumanPlaques #differential expression per HumanCluster

DiffExpHeat + DiffExpHeat_HumanPlaques

# Heatmap_with_genes_figure1Article ---------------------------------------

figure1genelist <- read_delim("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Figure1Genes.txt", delim = "\t")

All_Genes <- data.frame((rownames(Mouse_Obj.combined@assays$integrated@scale.data)))
colnames(All_Genes) <- "Gene"


Human_Top_Genes <- read_delim("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Human_Top_Genes_Cluster.csv", delim = ";") %>%
  separate(gene, c("gene", "Ensembl_gene"), sep = "-") %>%
  left_join(All_Genes, by = c("gene" = "Gene"), keep = TRUE) %>%
  drop_na(Gene) %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_logFC) %>%
  drop_na(Gene) %>%
  nest() 

for (i in 1:length(Human_Top_Genes$cluster)){
  print(i)
  Mouse_Obj.combined %<>%
    AddModuleScore(features = list(c(Human_Top_Genes$data[[i]]$Gene)), name = paste0("HumanPlaque", i), nbin = 24, ctrl = 100, assay = "integrated", slot = "scale.data", seed = 1)
  print(VlnPlot(Mouse_Obj.combined, features = paste("HumanPlaque", i ,"1",sep = "")))
}

Mouse_Obj.combined[["Module"]] <- CreateAssayObject(data = t(x = Seurat::FetchData(object = Mouse_Obj.combined, vars = c('HumanPlaque11', 'HumanPlaque21', 'HumanPlaque31', 'HumanPlaque41', 'HumanPlaque51')))) #adding module score to assay

DiffExpHeatFig1_data_mod <- DoHeatmap(Mouse_Obj.combined, features = c('HumanPlaque11', 'HumanPlaque21', 'HumanPlaque31', 'HumanPlaque41', 'HumanPlaque51' ) , assay = "Module", slot = "data")  + scale_fill_gradientn( colors = colorpanel(100,"grey100","grey90","darkgreen")) + labs(title = "Differential Gene Expression Mouse Dataset with genes from fig1 article")
DiffExpHeatFig1_data_mod

DiffExpHeatFig1_data_mod <- DoHeatmap(Mouse_Obj.combined, features = c('HumanPlaque11', 'HumanPlaque21', 'HumanPlaque31', 'HumanPlaque41', 'HumanPlaque51' ) , assay = "Module", slot = "data", group.by ="SeuratIntegrated")  + scale_fill_distiller(palette = "Greens", direction = 1) + labs(title = "Differential Gene Expression Mouse Dataset with genes from fig1 article")
DiffExpHeatFig1_data_mod

DiffExpHeatFig1_data <- DoHeatmap(Mouse_Obj.combined, features = Human_Top_Genes$gene)  + scale_fill_gradientn( colors = colorpanel(100,"grey100","grey90","darkgreen")) + labs(title = "Differential Gene Expression Mouse Dataset with genes from fig1 article")
DiffExpHeatFig1_data

DiffExpHeatFig1_data_HumanPlaq <- DoHeatmap(Mouse_Obj.combined, features = Human_Top_Genes$gene, group.by = "SeuratIntegrated")  + scale_fill_gradientn( colors = colorpanel(100,"grey100","grey90","darkgreen")) + labs(title = "Differential Gene Expression Mouse Dataset with genes from fig1 article")
DiffExpHeatFig1_data_HumanPlaq

DiffExpHeatFig1_data + DiffExpHeatFig1_data_HumanPlaq

DiffExpHeatFig1 <- DoHeatmap(Mouse_Obj.combined, features = figure1genelist$GenesSymbol) + scale_fill_gradientn( colors = colorpanel(100,"grey100","grey90","darkgreen"))+ labs(title = "Differential Gene Expression Mouse Dataset with genes from fig1 article")
DiffExpHeatFig1

DiffExpHeatFig1_HumanPlaq <- DoHeatmap(Mouse_Obj.combined, features = figure1genelist$GenesSymbol, group.by = "SeuratIntegrated")  + scale_fill_distiller(palette = "Greens", direction = 1) + labs(title = "Differential Gene Expression Mouse Dataset with genes from fig1 article")
DiffExpHeatFig1_HumanPlaq


# Clusterlist -------------------------------------------------------------

clusterlist <- data.frame(cbind(rownames(Mouse_Obj.combined@meta.data), as.character(Mouse_Obj.combined@meta.data$seurat_clusters)))
colnames(clusterlist) <- c("categories", "SeuratCluster")

JakeLusis %<>%
  left_join(clusterlist, by = c("categories"))
TableClus <- table(JakeLusis$SeuratCluster, JakeLusis$PlaqueCormat)


write.table(as.data.frame(JakeLusis) 
             #%>%
               #dplyr::select(-c(data,boxplot,annotateBox)))
             , 
            "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569_JakeLusis_Table_SeurataClustering.csv", 
            sep = ',', row.names = F, col.names = T, quote = F)

summary(as.factor(clusterlist$SeuratCluster))
summary(as.factor(JakeLusis$PlaqueCormat))

# Pathway_analysis --------------------------------------------------------
library(fastmap)
BiocManager::install("fastmap")
install.packages('fastmap', repos = c('http://rforge.net', 'http://cran.rstudio.org'),
                 type = 'source')
FindAllMarkers(Mouse_Obj.combined, only.pos= TRUE, min.pct = 0.1, logfc.threshold = 0.1) %>%
  group_by(cluster) %>%
  slice_max(n = 1500, order_by = avg_log2FC)  -> enrichment_paths 

Orthologs_Test <- getLDS( mart = Ens_Human,
                           attributes = 'hgnc_symbol',
                           martL = Ens_Mouse,
                           attributesL = c("entrezgene_id", "description"),
                           filters = 'hgnc_symbol',
                           values = enrichment_paths$gene) 

Orthologs_Test$NCBI.gene..formerly.Entrezgene..ID <- as.character(Orthologs_Test$NCBI.gene..formerly.Entrezgene..ID)

enrichment_paths %>%
  left_join(Orthologs_Test, by = c("gene"= "HGNC.symbol")) %>%
  group_by(cluster) %>%
  nest()%>%
  dplyr::mutate(pathway = map(data, ~enrichPathway(gene = .x$NCBI.gene..formerly.Entrezgene..ID, pvalueCutoff=0.05, organism = "mouse")),
                pathway_res = map(pathway, ~.x@result)) -> PathwayAnalysis

PathwayAnalysis %>% 
  dplyr::select(cluster, pathway_res) %>%
  unnest(pathway_res) %>%
  dplyr::mutate(geneID_lst = map(geneID, ~strsplit(.x, "/")))  %>%
  dplyr::select(ID, cluster, Description, GeneRatio, p.adjust, geneID_lst) %>%
  group_by(cluster) %>%
  slice_min(n = 15, order_by = p.adjust) %>%
  nest() -> Top_Pathways

Top_Pathways %>%
  unnest(cols = data) %>%
  group_by(cluster) %>%
  slice_min(n = 2, order_by = p.adjust) %>%
  nest() -> 
  

write.table(Top_n2_pathways %<>%
              unnest(cols = data) %>%
              ungroup() %>%
              dplyr::select(ID), 
            "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569_Topn2paths.csv", 
            sep = ',', row.names = F, col.names = T, quote = F)

for (cluster in 1:length(Top_Pathways$cluster)){
  print(cluster)
  for (pathways in 1:length(Top_Pathways$data[[cluster]]$geneID_lst)){
    print(Top_Pathways$data[[cluster]]$Description[[pathways]])
    Path_Genes <- Top_Pathways$data[[cluster]]$ID[[pathways]][[1]]
    
    Path_Genes_HGNC <- getLDS( mart = Ens_Mouse,
                               attributes = 'reactome',
                               martL = Ens_Human,
                               attributesL = c("hgnc_symbol"),
                               filters = 'reactome',
                               values = Path_Genes) 
    print(Path_Genes_HGNC$HGNC.symbol)
    if (length(Path_Genes_HGNC$HGNC.symbol > 1)){
      print("TRUE")
      Mouse_Obj.combined %<>%
        AddModuleScore(features = list(c(Path_Genes_HGNC$HGNC.symbol)), name = paste0( make.names(Top_Pathways$data[[cluster]]$Description[[pathways]]), ".cluster.", cluster), nbin = 24, ctrl = 100, assay = "integrated", slot = "scale.data", seed = 1)
      print(VlnPlot(Mouse_Obj.combined, features = paste0(make.names(Top_Pathways$data[[cluster]]$Description[[pathways]]), ".cluster.", cluster, 1))) } 
  }
} 


Path_select <- data.frame(Id =  c("R-MMU-9020591", "R-MMU-446353", "R-MMU-203615", "R-MMU-165159", "R-MMU-917937", "R-MMU-71240", "R-MMU-168256", "R-MMU-168898", "R-MMU-913531", "R-MMU-168249", "R-MMU-6798695", "R-MMU-8983432", "R-MMU-70171")) #some selected pathways

for (j in 1: length(Path_select$Id)){
  Path_Name <- query(Path_select[j,])$displayName
  print(Path_Name)
 print(Path_select[j,])
    pathway=c(Path_select[j,])
  path_genes <- getLDS( mart = Ens_Mouse,
                        attributes = 'reactome',
                        martL = Ens_Human,
                        attributesL = c("hgnc_symbol"),
                        filters = 'reactome',
                        values = pathway)
  print(path_genes$HGNC.symbol)
     Mouse_Obj.combined %<>%
       AddModuleScore(features = list(c(path_genes$HGNC.symbol)) , name = paste0( make.names(Path_Name), "."), nbin = 24, ctrl = 100, assay = "integrated", slot = "scale.data", seed = 1)
     print(VlnPlot(Mouse_Obj.combined, features =paste0( make.names(Path_Name),".",  "1")))
}


Path_Names <- colnames(Mouse_Obj.combined@meta.data) %>%
  str_extract_all(pattern = ".+cluster\\.\\d1") %>%
  unlist() #extract all the pathways from metadata

Mouse_Obj.combined[["Module"]] <- CreateAssayObject(data = t(x = Seurat::FetchData(object = Mouse_Obj.combined, vars = Path_Names))) 

DiffExpHeatFig1_data_mod <- DoHeatmap(Mouse_Obj.combined, features = c(Path_Names) , assay = "Module", slot = "data")  + scale_fill_gradient2(low = "black", mid = "grey90", high = "darkgreen", midpoint = 0) + labs(title = "Differential Gene Expression Mouse Dataset with genes from fig1 article")
DiffExpHeatFig1_data_mod

DiffExpHeatFig1_data_mod_HumanPlaq <- DoHeatmap(Mouse_Obj.combined, features = c(Path_Names) , assay = "Module", slot = "data", group.by = "SeuratIntegrated")  + scale_fill_gradient2(low = "black", mid = "grey90", high = "darkgreen", midpoint = 0) + labs(title = "Differential Gene Expression Mouse Dataset with genes from fig1 article")
DiffExpHeatFig1_data_mod_HumanPlaq

Path_Names %<>%
  as.data.frame() %>%
  dplyr::mutate(VlnPlot = map(Path_Names, ~VlnPlot(object = Mouse_Obj.combined, features = .x)))


Vln_plt_clus1 <- Path_Names$VlnPlot[[1]] + labs(title = "Fatty acid metabolism Cluster A") + stat_compare_means( method = "anova", label.y.npc = "top", label.x.npc = "middle") + stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  theme(plot.title = element_text(face = "plain"))
Vln_plt_clus2 <- Path_Names$VlnPlot[[16]] + labs(title = "Extracellular matrix organization Cluster B") + stat_compare_means( method = "anova", label.y.npc = "top", label.x.npc = "middle") + stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  theme(plot.title = element_text(face = "plain"))
Vln_plt_clus3 <- Path_Names$VlnPlot[[32]] + labs(title = "Muscle contraction Cluster C") + stat_compare_means( method = "anova", label.y.npc = "top", label.x.npc = "middle") + stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  theme(plot.title = element_text(face = "plain"))
Vln_plt_clus4 <- Path_Names$VlnPlot[[50]] + labs(title = "Neuronal Pathway Cluster D") + stat_compare_means( method = "anova", label.y.npc = "top", label.x.npc = "middle") + stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  theme(plot.title = element_text(face = "plain"))

Vln_Pathway_All_Clus <- grid.arrange(Vln_plt_clus1, Vln_plt_clus2, Vln_plt_clus3, Vln_plt_clus4, ncol = 4)

d1 <- dotplot(PathwayAnalysis$pathway[[1]], showCategory=10) + labs(title = "Pathways Cluster A")  + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))
d2 <- dotplot(PathwayAnalysis$pathway[[2]], showCategory=10) + labs(title = "Pathways Cluster B") + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))
d3 <- dotplot(PathwayAnalysis$pathway[[3]], showCategory=10) + labs(title = "Pathways Cluster C") + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))
d4 <- dotplot(PathwayAnalysis$pathway[[4]], showCategory=10) + labs(title = "Pathways Cluster D") + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))

d1 + d2 + d3 + d4

d1_poster <- dotplot(PathwayAnalysis$pathway[[1]], showCategory=4) + labs(title = "Pathways Cluster A")  + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))
d2_poster <- dotplot(PathwayAnalysis$pathway[[2]], showCategory=4) + labs(title = "Pathways Cluster B") + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))
d3_poster <- dotplot(PathwayAnalysis$pathway[[3]], showCategory=4) + labs(title = "Pathways Cluster C") + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))
d4_poster <- dotplot(PathwayAnalysis$pathway[[4]], showCategory=4) + labs(title = "Pathways Cluster D") + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))

d1_poster + d2_poster + d3_poster + d4_poster
Pathway_All_Clus <- grid.arrange(d1, d2, d3, d4)
Pathway_All_Clus
library(gridExtra)
Pathway_All_Clus_Poster <- grid.arrange(d1_poster , d2_poster , d3_poster , d4_poster)
#d1_poster <- dotplot(PathwayAnalysis$pathway[[1]][1:4,], showCategory=10) + labs(title = "Pathways Cluster A")  + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))
#d2_poster <- dotplot(PathwayAnalysis$pathway[[2]], showCategory=10) + labs(title = "Pathways Cluster B") + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))
#d3_poster <- dotplot(PathwayAnalysis$pathway[[3]], showCategory=10) + labs(title = "Pathways Cluster C") + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))
#d4_poster <- dotplot(PathwayAnalysis$pathway[[4]], showCategory=10) + labs(title = "Pathways Cluster D") + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))

#d1 + d2 + d3 + d4
Path_Names_Select <- colnames(Mouse_Obj.combined@meta.data) %>%
  str_extract_all(pattern = ".+\\.1$") %>%
  unlist() %>%
  as.data.frame() 
Path_Names_Select %<>%
  dplyr::mutate(VlnPlot = map(Path_Names_Select[,1], ~VlnPlot(object = Mouse_Obj.combined, features = .x) + stat_compare_means( method = "anova", label.y.npc = "top", label.x.npc = "middle") + stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") ),
                VlnPlot_HumanPlq = map(Path_Names_Select[,1], ~VlnPlot(object = Mouse_Obj.combined, features = .x, group.by = "SeuratIntegrated")+ stat_compare_means( method = "anova", label.y.npc = "top", label.x.npc = "middle") + stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") ))


Vln_select1 <- Path_Names_Select$VlnPlot[[1]] + labs(title = "Interleukin 12 signalling") +
  theme(plot.title = element_text(face = "plain"))
Vln_select2 <- Path_Names_Select$VlnPlot[[7]] + labs(title = "Immune system") +
  theme(plot.title = element_text(face = "plain"))
Vln_select3 <- Path_Names_Select$VlnPlot[[13]] + labs(title = "Glycolysis") +
  theme(plot.title = element_text(face = "plain"))
Vln_select4 <- Path_Names_Select$VlnPlot[[4]] + labs(title = "MTOR signalling") +
  theme(plot.title = element_text(face = "plain"))
Vln_select5 <- Path_Names_Select$VlnPlot[[11]] + labs(title = "Neutrophil degranulation") +
  theme(plot.title = element_text(face = "plain"))
Vln_select6 <- Path_Names_Select$VlnPlot[[5]] + labs(title = "Iron uptake and transport") +
  theme(plot.title = element_text(face = "plain"))
Vln_select7 <- Path_Names_Select$VlnPlot[[2]]  +
  theme(plot.title = element_text(face = "plain"))



Vln_select1 + Vln_select2 + Vln_select3 + Vln_select4 + Vln_select5 + Vln_select6

Vln_Pathway_Select_Clus <- grid.arrange(Vln_select1, Vln_select2, Vln_select3, Vln_select4,Vln_select5,  Vln_select6, ncol = 3)
Vln_Pathway_Select_Clus


Heatmap_Paths <- Top_Pathways[c(1,16,32, 50),]
Top_Pathways %>%
  unnest (cols = c(data)) -> Heatmap_Paths
Heatmap_Paths[c(1,16,32, 50),] %>%
  dplyr::mutate(Heatmap_genes= map(ID, ~ DoHeatmap(Mouse_Obj.combined, features = getLDS( mart = Ens_Mouse,
                                                                                 attributes = 'reactome',
                                                                                 martL = Ens_Human,
                                                                                 attributesL = c("hgnc_symbol"),
                                                                                 filters = 'reactome',
                                                                                values = .x)$HGNC.symbol) + scale_fill_gradientn( colors = colorpanel(100,"grey100","grey90","darkgreen")))) -> a
a$Heatmap_genes


# Violin_FeaturePlots -----------------------------------------------------

LPIN <- FeaturePlot(Mouse_Obj.combined, c("LPIN1")) + ggtitle("LPIN1 Cluster #0")
C1QA <- FeaturePlot(Mouse_Obj.combined, c("C1QA")) + ggtitle("C1QA Cluster #1")
ACTA2 <- FeaturePlot(Mouse_Obj.combined, c("ACTA2")) + ggtitle("ACTA2 Cluster #2")
CD74 <- FeaturePlot(Mouse_Obj.combined, c("NOS1")) + ggtitle("NOS1 Cluster #3")
MYH11 <- FeaturePlot(Mouse_Obj.combined, c("MYH11")) + ggtitle("MYH11 Cluster #4")

LPIN_vln <- VlnPlot(Mouse_Obj.combined, c("LPIN1")) + ggtitle("LPIN1 Cluster #0")
C1QA_vln <- VlnPlot(Mouse_Obj.combined, c("C1QA")) + ggtitle("C1QA Cluster #1")
ACTA2_vln <- VlnPlot(Mouse_Obj.combined, c("ACTA2")) + ggtitle("ACTA2 Cluster #2")
CD74_vln <- VlnPlot(Mouse_Obj.combined, c("NOS1")) + ggtitle("NOS1 Cluster #3")
MYH11_vln <- VlnPlot(Mouse_Obj.combined, c("MYH11")) + ggtitle("MYH11 Cluster #4")

HumanViolinplots <- grid.arrange(LPIN_vln , C1QA_vln , ACTA2_vln , CD74_vln , MYH11_vln, ncol = 2, top = paste("Human Genes violin plots"))
HumanFeatures <- grid.arrange(LPIN , C1QA , ACTA2 , CD74 , MYH11, ncol = 2, top = paste("Human Genes Feature plots"))

LIPA <- FeaturePlot(Mouse_Obj.combined, c("LIPA"),pt.size = 1.5) + ggtitle("LIPA Cluster #0")
SERPINE2 <- FeaturePlot(Mouse_Obj.combined, c("SERPINE2"),pt.size = 1.5) + ggtitle("SERPINE2 Cluster #2")
ACTA1 <- FeaturePlot(Mouse_Obj.combined, c("ACTA1"),pt.size = 1.5) + ggtitle("ACTA1 Cluster #2")
NPY <- FeaturePlot(Mouse_Obj.combined, c("NPY"),pt.size = 1.5) + ggtitle("NPY Cluster #3")

LIPA_vln <- VlnPlot(Mouse_Obj.combined, c("LIPA"),pt.size = 1.5) + ggtitle("LIPA Cluster #0")
SERPINE2_vln <- VlnPlot(Mouse_Obj.combined, c("SERPINE2"),pt.size = 1.5) + ggtitle("SERPINE2 Cluster #2")
ACTA1_vln <- VlnPlot(Mouse_Obj.combined, c("ACTA1"),pt.size = 1.5) + ggtitle("ACTA1 Cluster #2")
NPY_vln <- VlnPlot(Mouse_Obj.combined, c("NPY"),pt.size = 1.5) + ggtitle("NPY Cluster #3")

MouseViolinplots <- arrangeGrob(LIPA_vln, SERPINE2_vln , ACTA1_vln , NPY_vln, top = textGrob("Mouse Genes violin Plots", gp = gpar(fontsize = 17, fontface = 'bold')), ncol = 1, nrow = 4)
MouseFeatures <- grid.arrange(LIPA, SERPINE2 , ACTA1 , NPY, top = textGrob("Mouse Genes Feature Plots", gp = gpar(fontsize = 17, fontface = 'bold')), ncol = 1, nrow = 4)

# Frequent_Strain ---------------------------------------------------------
main_strain = c(  "BXD/RwwJ", 
                  "BXD/TyJ",
                  "AXB/PgnJ", 
                  "BXA/PgnJ",
                  "CXB/ByJ",
                  "BXH/TyJ")
JakeLusis %>%
  dplyr::select(c(Strain, SeuratCluster, categories,PlaqueCormat )) %>%
  dplyr::mutate(cleaned_strains = str_replace_all(Strain, c( "BXD\\d+\\/RwwJ" = "BXD/RwwJ", 
                                                             "^(BXD\\d+\\/TyJ)" = "BXD/TyJ",
                                                             "AXB(\\d+.|\\d+)\\/PgnJ" = "AXB/PgnJ", 
                                                             "BXA\\d+\\/PgnJ" = "BXA/PgnJ",
                                                             "CXB\\d+\\/ByJ" = "CXB/ByJ",
                                                             "BXH\\d+\\/TyJ" = "BXH/TyJ"))) %>%
  dplyr::filter(str_detect(cleaned_strains, paste(main_strain, collapse = "|"))) -> Selected_Strains

JakeLusis %>%
  dplyr::select(c(Strain, SeuratCluster, categories, PlaqueCormat)) %>%
  mutate_at('SeuratCluster', as.character) %>%
  group_by(Strain)  %>%
  nest() %>%
  dplyr::mutate(cluster1 = map(data, ~sum(.x$SeuratCluster == "1")),
                cluster2 = map(data, ~sum(.x$SeuratCluster == "2")),
                cluster3 = map(data, ~sum(.x$SeuratCluster == "3")),
                cluster4 = map(data, ~sum(.x$SeuratCluster == "4")))  %>%
  dplyr::select(-c(data)) -> Strain_per_cluster

  
Selected_Strains %>%
  dplyr::select(!PlaqueCormat) %>%
  group_by(SeuratCluster) %>%
  nest() %>%
  dplyr::mutate(count = map(data, ~as.vector(.x$cleaned_strains))) %>%
  dplyr::mutate(!!paste0(main_strain[1]) := str_count(count, paste(main_strain[1])),
                !!paste0(main_strain[2]) := str_count(count, paste(main_strain[2])),
                !!paste0(main_strain[3]) := str_count(count, paste(main_strain[3])),
                !!paste0(main_strain[4]) := str_count(count, paste(main_strain[4])),
                !!paste0(main_strain[5]) := str_count(count, paste(main_strain[5])),
                !!paste0(main_strain[6]) := str_count(count, paste(main_strain[6]))) -> Main_Strain_Count_SeuratClus
Selected_Strains %>%
  dplyr::select(!SeuratCluster) %>%
  group_by(PlaqueCormat) %>%
  nest() %>%
  dplyr::mutate(count = map(data, ~as.vector(.x$cleaned_strains))) %>%
  dplyr::mutate(!!paste0(main_strain[1]) := str_count(count, paste(main_strain[1])),
                !!paste0(main_strain[2]) := str_count(count, paste(main_strain[2])),
                !!paste0(main_strain[3]) := str_count(count, paste(main_strain[3])),
                !!paste0(main_strain[4]) := str_count(count, paste(main_strain[4])),
                !!paste0(main_strain[5]) := str_count(count, paste(main_strain[5])),
                !!paste0(main_strain[6]) := str_count(count, paste(main_strain[6]))) -> Main_Strain_Count_PlaqueCor

Selected_Strains %>%
  select(cleaned_strains, categories) -> Main_Strain_list

Mouse_Obj.combined@meta.data %<>%
  rownames_to_column(var = "categories") %>% 
  left_join(Main_Strain_list,  by = c("categories")) %>%
  column_to_rownames(var = "categories")

MouseUmapSeuratSelectedStrain <- DimPlot(Mouse_Obj.combined, reduction = "umap", group.by = "cleaned_strains", pt.size = 1.5)
MouseUmapSeuratSelectedStrain

MouseUmap + MouseUmapSeuratSelectedStrain
JakeLusis %>%
  group_by(Strain) %>%
  nest() %>%
  dplyr::mutate(seuratclus = map(data, ~.x %>% 
                                   group_by(SeuratCluster) %>% 
                                  dplyr::summarize(n = n()) %>% 
  top_n(1))) %>%
  dplyr::mutate(seuratclus_Mean = map(data, ~weighted.mean(as.numeric(.x$SeuratCluster)))) -> SummarizedJakeLusis


ggsave(paste0(DataName, "JakelusisMeasurements_SeuratClusters.pdf"), Measurement_plots, width = 100, height = 80, units = "cm")

Idents(Mouse_Obj.combined) <- Mouse_Obj.combined@meta.data$SeuratIntegrated

FindAllMarkers(Mouse_Obj.combined, only.pos= TRUE, min.pct = 0.1, logfc.threshold = 0.1) %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC) -> conserved_markers



Orthologs_Conserved <- getLDS( mart = Ens_Human,
                          attributes = 'hgnc_symbol',
                          martL = Ens_Mouse,
                          attributesL = c("entrezgene_id", "description"),
                          filters = 'hgnc_symbol',
                          values = conserved_markers$gene) 

Orthologs_Conserved$NCBI.gene..formerly.Entrezgene..ID <- as.character(Orthologs_Conserved$NCBI.gene..formerly.Entrezgene..ID)

conserved_markers %>%
  left_join(Orthologs_Conserved, by = c("gene"= "HGNC.symbol")) %>%
  group_by(cluster) %>%
  nest()%>%
  dplyr::mutate(pathway = map(data, ~enrichPathway(gene = .x$NCBI.gene..formerly.Entrezgene..ID, pvalueCutoff=0.05, organism = "mouse")),
                pathway_res = map(pathway, ~.x@result)) -> PathwayAnalysis_conserved
library(gplots)
d1 <- dotplot(PathwayAnalysis_conserved$pathway[[1]], showCategory=10) + labs(title = "Pathways Cluster 4")  + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))
d2 <- dotplot(PathwayAnalysis_conserved$pathway[[2]], showCategory=10) + labs(title = "Pathways Cluster 1") + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))
d3 <- dotplot(PathwayAnalysis_conserved$pathway[[3]], showCategory=10) + labs(title = "Pathways Cluster 3") + scale_colour_gradientn( colors = colorpanel(100,"darkgreen","grey70"))

d1 + d2 + d3 

DiffExpHeat <- DoHeatmap(Mouse_Obj.combined, features = conserved_markers$gene)  + scale_fill_distiller(palette = "Greens", direction = 1) + labs(title = "Differential Gene Expression Mouse Dataset, MouseClusters")
Human_Mice_Aorta_Coronary_Obj@meta.data


# Markers -----------------------------------------------------------------
NCount_RNA_MouseAPOE <- VlnPlot(Human_Mice_Aorta_Coronary_Obj, features = "nCount_RNA", group.by = "cluster")

Integrated_carotid %>%
  left_join(Integrated_coronary) %>%
  distinct(gene_id , .keep_all = TRUE) %>%
  left_join(Healthy_Aorta, by = c("gene_id")) %>%
  distinct(gene_id , .keep_all = TRUE) %>%
  left_join(Integrated_Mouse_all_genes %>%
              rownames_to_column(var = "gene_id"), by = c("gene_id")) %>%
  distinct(gene_id , .keep_all = TRUE) %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'gene_id') %>%
  normalizeQuantiles() -> Human_Mice_Aorta_Coronary_dataframe

Human_Mice_Aorta_Coronary_dataframe %>%
  dplyr::mutate(across(everything(), rank)) %>%
  dplyr::mutate(across(everything(), function(x) 100 * x / max(x)))-> Ranked_Df_perc

Ranked_perc_Obj <- CreateSeuratObject(counts = Ranked_Df_perc)

Ranked_perc_Obj@meta.data %<>%
  rownames_to_column("Name") %>%
  dplyr::mutate(dataset = ifelse(grepl("ae",Name), paste0("Human_Caro"), 
                                 ifelse(grepl("GSM",Name), paste0("Mouse"),
                                        ifelse(grepl("GTEX",Name), paste0("Human_Aorta"),
                                               ifelse(grepl("UVA",Name), paste0("Human_Cor"), "-"))))) %>%
  left_join(Human_Mouse_Cats, by = c("Name" = "categories")) %>%
  left_join(SeuratClusters, by = c("Name" = "study_number")) %>%
  dplyr::mutate(cluster = coalesce(cluster, cats)) %>%
  mutate_at('cluster',  ~replace_na(., "Human_Aorta"))  %>%
  column_to_rownames("Name")


VEGFD <- VlnPlot(Human_Mice_Aorta_Coronary_Obj, features = "VEGFA", group.by = "cluster", slot = "scale.data") + NoLegend()
VCAM <- VlnPlot(Human_Mice_Aorta_Coronary_Obj, features = "VCAM1", group.by = "cluster", slot = "scale.data") + NoLegend()
ICAM <- VlnPlot(Human_Mice_Aorta_Coronary_Obj, features = "ICAM1", group.by = "cluster", slot = "scale.data") + NoLegend()
TGF <- VlnPlot(Human_Mice_Aorta_Coronary_Obj, features = "TGFB1", group.by = "cluster", slot = "scale.data") + NoLegend()
IL1B <- VlnPlot(Human_Mice_Aorta_Coronary_Obj, features = "IL1B", group.by = "cluster", slot = "scale.data") + NoLegend()
IL6 <- VlnPlot(Human_Mice_Aorta_Coronary_Obj, features = "IL6", group.by = "cluster", slot = "scale.data") + NoLegend()

LPIN1 <- VlnPlot(Ranked_perc_Obj, features = "LPIN1", group.by = "cluster", slot = "data") + NoLegend() + ylab("Percentile ranked\n gene expression") 
KYNU <- VlnPlot(Ranked_perc_Obj, features = "KYNU", group.by = "cluster", slot = "data") + NoLegend() + ylab("Percentile ranked\n gene expression")
CD14 <- VlnPlot(Ranked_perc_Obj, features = "CD14", group.by = "cluster", slot = "data") + NoLegend() + ylab("Percentile ranked\n gene expression")
MYH11 <- VlnPlot(Ranked_perc_Obj, features = "MYH11", group.by = "cluster", slot = "data") + NoLegend() + ylab("Percentile ranked\n gene expression")
ACTA2 <- VlnPlot(Ranked_perc_Obj, features = "ACTA2", group.by = "cluster", slot = "data") + NoLegend() + ylab("Percentile ranked\n gene expression")
NOS1 <- VlnPlot(Ranked_perc_Obj, features = "NOS1", group.by = "cluster", slot = "data") + NoLegend() + ylab("Percentile ranked\n gene expression")
SOD2 <- VlnPlot(Ranked_perc_Obj, features = "SOD2", group.by = "cluster", slot = "data") + NoLegend() + ylab("Percentile ranked\n gene expression")

LPIN1 + KYNU + CD14 + MYH11 + ACTA2 + NOS1 + SOD2

VEGFA + VEGFD

VEGFA + VCAM + ICAM + TGF + IL1B + IL6
Markers_atherosclerosis <- grid.arrange(VEGFD, VCAM , ICAM , IL1B, IL6, top = textGrob("Important Markers Atherosclerosis field", gp = gpar(fontsize = 17, fontface = 'bold')), ncol = 2)

Markers_atherosclerosis2 <- grid.arrange(LPIN1, KYNU , CD14 , MYH11, ACTA2, NOS1, SOD2,  top = textGrob("Important Markers Atherosclerosis field", gp = gpar(fontsize = 17, fontface = 'bold')), ncol = 2)
Markers_atherosclerosis2
library(gridExtra)



markers <- FindMarkers(Human_Mice_Aorta_Coronary_Obj, indent.1 = "Carotid_1",   test.use = "roc")

# Measurement boxplots ----------------------------------------------------
#making a list with possible measurement traits
Measurements <- read_delim("C:/Users/brouw/Downloads/pgen.1005711.s014.csv") %>%
  group_by(trait_name) %>%
  nest() %>%
  dplyr::select(trait_name)

write.table(Measurements, "C:/Users/brouw/Downloads/Measurements.txt", sep = ",", row.names = F)

#pdf with boxplts
Measurements_JakeLusis <- read_delim("C:/Users/brouw/Downloads/pgen.1005711.s014_Num.csv") %>%
  dplyr::select(-c(...7)) %>%
  group_by( Maternal_strain ) %>%
  nest() %>% 
  unnest(cols = data) 
my_comparisons = list(c("1", "2"), c("2", "3"), c("3", "4"), c("1", "4"))

JakeLusis %>%
  dplyr::select(c(categories, Strain, SeuratCluster)) %>%
  left_join(Measurements_JakeLusis, by = c("Strain" = "Maternal_strain")) %>% 
  group_by(trait_name) %>%
  nest() %>%
  dplyr::mutate(anova = purrr::map(data, ~summary(aov(.x$SeuratCluster ~.x$avg) ))) %>%
  #dplyr::mutate(corrected_pval = purrr::map(anova, ~p.adjust(.x$"Pr(>F)", "fdr"))) %>%
  dplyr::mutate(boxplot = purrr::map2(data, trait_name, ~ ggplot(data = .x , aes(x = SeuratCluster, y = avg, fill = factor(SeuratCluster))) + 
                                        geom_boxplot() +
                                        scale_fill_brewer(palette = "Paired") + 
                                        theme(legend.position = "none") +
                                        stat_compare_means( method = "anova", label.y.npc = "top", label.x.npc = "middle") +
                                        stat_compare_means(comparisons = my_comparisons) +
                                        xlab(.y) )) -> Measurement_plots


Measurement_plots <- marrangeGrob(Measurement_plots$boxplot,  nrow = 4, ncol = 5 )

Plots_top_Measurements <- grid.arrange(Measurement_plots$boxplot[1][[1]], Measurement_plots$boxplot[14][[1]], Measurement_plots$boxplot[78][[1]], Measurement_plots$boxplot[76][[1]], Measurement_plots$boxplot[21][[1]])

#table with measurements
JakeLusis %>%
  dplyr::select(c(categories, Strain, SeuratCluster)) %>%
  left_join(Measurements_JakeLusis, by = c("Strain" = "Maternal_strain")) %>%
  gtsummary::select(categories, trait_name, SeuratCluster, avg) %>%
  group_by(categories) %>%
  spread(trait_name, avg) %>%
  column_to_rownames(var = "categories") %>%
  tbl_summary(by = SeuratCluster,
              missing = "no") %>%
  add_overall(last = "true") %>%
  modify_header(label ~ "**Trait measured**") %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4") ~ "**Mouse Clusters**") %>%
  add_p(test = everything() ~"aov")  -> trait_table

trait_table %>% 
  as_hux_xlsx("Measurement_table.xlsx")

