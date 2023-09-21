#Human Carotids Human Coronary integrated for pathway correlation

#end product
Integrated_carotid <- read_delim("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Integrated_coronary_carotids_carotidsamples.txt", show_col_types = FALSE)
Integrated_coronary <- read_delim("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Integrated_coronary_carotids_coronarysamples.txt", show_col_types = FALSE)


Human_Car_Pth <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Human/Human1/bulk_RNAseq_raw_counts.txt.minRib.txt.PC.txt"
Human_Cor_Pth <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Human/Human2/genecountsraw_no_aorta.txt.PC.txt.minRib2.txt"
Coronary_Clus_Pth <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Seurat_clusters_Michal_v13.txt.parsed.merged_clusters_clean_Human2Cluster1.txt"


Human_Cor <- read_delim(Human_Cor_Pth, show_col_types = FALSE) %>%
  dplyr::select(-c(ENSID)) %>%
  dplyr::rename(gene_id =  Description) %>%
  remove_rownames() %>%
  group_by(gene_id) %>%
  dplyr::summarize_if(is.numeric, mean) %>% #every duplicate will be summarized
  drop_na() %>%
  column_to_rownames(var = 'gene_id') %>% #no further normalization, data is already correct form
  normalizeQuantiles()  %>% 
  rownames_to_column(var = "gene_id") 

Human_Caro <- read_delim(Human_Car_Pth, show_col_types = FALSE) %>%
  dplyr::select(-c(start,end,ensembl_gene_id,  strand)) %>%
  dplyr::rename(gene_id = symbol) %>%
  remove_rownames() %>%
  group_by(gene_id) %>%
  dplyr::summarize_if(is.numeric, mean) %>% 
  drop_na() %>%
  column_to_rownames(var = 'gene_id') %>%
  mutate_all( ~ifelse(.x > 4095, 4095, .x)) %>%
  mutate_all( ~-4096*log(1-(.x/4096))) %>% #data normalization
  normalizeQuantiles()  %>%
  rownames_to_column(var = "gene_id")

Coronary_Clus <- read_delim(Coronary_Clus_Pth) 


#adding coronary clustering

Coronary_cluster0 <- Coronary_Clus[which(Coronary_Clus$cluster == 0), ] #list of classes
Coronary_cluster1 <- Coronary_Clus[which(Coronary_Clus$cluster == 1), ]
Coronary_cluster2 <- Coronary_Clus[which(Coronary_Clus$cluster == 2), ]
Coronary_cluster3 <- Coronary_Clus[which(Coronary_Clus$cluster == 3), ]
Coronary_cluster4 <- Coronary_Clus[which(Coronary_Clus$cluster == 4), ]

#joining the two matrixes
Human_both <- Human_Cor %>%
  left_join(Human_Caro, by = c("gene_id")) %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = "gene_id") %>%
  normalizeQuantiles()

Human_Obj <- CreateSeuratObject(counts = Human_both)

Human_Obj@meta.data %<>%
  rownames_to_column("Names") %>%
  dplyr::mutate(dataset = ifelse(grepl("ae",variable), paste0("Human_Caro"), 
                                               ifelse(grepl("UVA",variable), paste0("Human_Cor"), "-"))) %>%
  dplyr::mutate(Plaque_Cat = ifelse(Names %in% cluster0$study_number, paste0("0"),
                                    ifelse(Names %in% cluster1$study_number,paste0("1"),
                                           ifelse(Names %in% cluster2$study_number, paste0("2"),
                                                  ifelse(Names %in% cluster3$study_number, paste0("3"),
                                                         ifelse(Names %in% cluster4$study_number, paste0("4"), 
                                                                ifelse(Names %in% Coronary_cluster0$study_number, paste0("0"),
                                                                      ifelse(Names %in% Coronary_cluster1$study_number, paste0("1"),
                                                                              ifelse(Names %in% Coronary_cluster2$study_number, paste0("2"),
                                                                                     ifelse(Names %in% Coronary_cluster3$study_number, paste0("3"), 
                                                                                            ifelse(Names %in% Coronary_cluster4$study_number, paste0("4"), NA))))))))))) %>%
  column_to_rownames("Names")

Human_Obj.list <- SplitObject(Human_Obj, split.by = "dataset")

Human_Obj.list <- lapply(X = Human_Obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000000)
})

features <- SelectIntegrationFeatures(object.list = Human_Obj.list)

Human_Obj.anchors <- FindIntegrationAnchors(object.list = Human_Obj.list, anchor.features = features)
Human_Obj.combined <- IntegrateData(anchorset = Human_Obj.anchors)

DefaultAssay(Human_Obj.combined) <- "integrated"

#add the categories column
Human_Obj.combined@meta.data %<>%
  rownames_to_column(var = 'Name') %>%
  left_join(Human_Mouse_Cats, by = c("Name" = "categories")) %>%
  left_join(SeuratClusters, by = c("Name" = "study_number")) %>%
  dplyr::mutate(cluster = coalesce(cluster, cats)) %>%
  mutate_at('cluster',  ~replace_na(., "Human_Aorta"))  %>%
  column_to_rownames("Name")

Human_Obj.combined %<>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>% 
  RunUMAP( dims = 1:12) 

Human_Obj.combined.list <- SplitObject(Human_Obj.combined, split.by = "dataset")

as.data.frame(GetAssayData(object = Human_Obj.combined.list[[1]], assay = "RNA", slot = "data")) -> Integrated_coronary
write.table(as.data.frame(Integrated_coronary) %>%
            rownames_to_column("gene_id")
            , 
            "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Integrated_coronary_carotids_coronarysamples.txt", 
            sep = ',', row.names = F, col.names = T, quote = F)

as.data.frame(GetAssayData(object = Human_Obj.combined.list[[2]], assay = "RNA", slot = "data")) -> Integrated_carotid

write.table(as.data.frame(Integrated_carotid) %>%
              rownames_to_column("gene_id")
            , 
            "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Integrated_coronary_carotids_carotidsamples.txt", 
            sep = ',', row.names = F, col.names = T, quote = F)
# clustering and correlation ----------------------------------------------
p3 <- DimPlot(Human_Obj.combined, reduction = "umap",split.by = c('dataset') , group.by = c( "Plaque_Cat") , pt.size = 1.5)
p3

cormat <- as.data.frame(stats::cor(as.matrix(Seurat::GetAssayData(Human_Obj.combined, assay = "integrated", slot = "data")))) %>%
  tibble::rownames_to_column(var = "Names") 

Human <- tibble(.rows = length(cormat)-1)

Human %<>%
  dplyr::mutate(Names = cormat$Names) %>%
  left_join(cormat, by = c("Names")) %>%
  dplyr::mutate(
    Plaque_Cat = ifelse(Names %in% cluster0$study_number, 0, 
                        ifelse(Names %in% cluster1$study_number,1, 
                               ifelse(Names %in% cluster2$study_number, 2, 
                                      ifelse(Names %in% cluster3$study_number, 3, 
                                             ifelse(Names %in% cluster4$study_number, 4, NA)))))) %>% #changed the NA from Dataname
  na.omit() 

Human %>%
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
                "4" = map_dbl(cor_plaque4, ~median(.x$Correlation))) %>%
  dplyr::select(-c(cor_plaque0, cor_plaque1, cor_plaque2, cor_plaque3, cor_plaque4)) %>%
  pivot_longer(-c(categories, data), names_to = "PlaqueCormat") %>%
  group_by(categories) %>% 
  slice_max(value) %>% 
  ungroup() %>%
  left_join(SeuratClusters, by = c( "categories" = "study_number" )) %>%
  dplyr::mutate(comparison = ifelse(PlaqueCormat > cluster, "X",
                                    ifelse(PlaqueCormat < cluster, "X", ""))) %>%
  dplyr::select(-c(data)) -> Mice_CormatPlaque

Human_plots <- Human %>%
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

Tissue_compare <- Human_plots[657:818,] 
Tissue_compare

Human_Obj.combined@meta.data %<>%
  rownames_to_column(var = "categories") %>%
  left_join(Tissue_compare, by = "categories") %>%
  dplyr::select(-c(data,boxplot,value,cluster,comparison,annotateBox)) %>%
  column_to_rownames(var = 'categories')

p3 <- DimPlot(Human_Obj.combined, reduction = "umap",split.by = c('dataset') , group.by = c( "PlaqueCormat") , pt.size = 1.5)
p3

tableJake <- table(Tissue_compare$PlaqueCormat)
tableJake

