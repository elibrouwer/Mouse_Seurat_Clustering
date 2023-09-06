#Script to remove the APOE Mouse Batch effect
pacman::p_load("GEOquery", "Seurat")

# Tables that are produced by this script ---------------------------------
Strain_Batch <- read_csv("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569_Strain_Batch.csv")

Integrated_Mouse <- read_csv("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569_Corrected_Batch.csv") %>%
  column_to_rownames(var = "gene_id")
Integrated_Mouse_all_genes <- read_csv("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569_Corrected_Batch_All_Genes.csv") %>%
  column_to_rownames(var = "gene_id")

# Remove Batch Effect -----------------------------------------------------

## making a table to add strains and batch number
gse=getGEO(filename="C:/Users/brouw/Downloads/GSE66569_series_matrix (2).txt.gz")
as.data.frame(str_split(gse@phenoData@data$supplementary_file, pattern = "/", simplify = TRUE)) %>%
  dplyr::select(V9) -> col_list
col_list$Name <- col_list$V9
col_list$Strain <- gse@phenoData@data$`strain:ch1`

as.data.frame(str_split(col_list$Name, pattern = "_", simplify = TRUE)) %>%
  group_by(V2) %>%
  nest() %>%
  rownames_to_column(var = "Batch") %>%
  unnest(cols = c(data)) %>%
  dplyr::mutate(id = paste0(V1, "_", V2, "_", V3)) -> c
col_list %>%
  left_join(c, by = c("Name" = "id")) %>%
  dplyr::select(c(Name, Strain, Batch)) -> Strain_Batch

table(factor(Strain_Batch$Batch))

#Strain_Batch --> table with the different batches
write.table((as.data.frame(Strain_Batch)), 
            "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569_Strain_Batch.csv", 
            sep = ',', row.names = F, col.names = T, quote = F)

# correction for batches --------------------------------------------------
APOE_Mouse <- read_delim(MousePath) #RMA normalized (log2)

orthologs_mouse <- getLDS( mart = Ens_Mouse,
                           attributes = 'ensembl_gene_id',
                           martL = Ens_Human,
                           attributesL = c( "hgnc_symbol"),
                           filters = "ensembl_gene_id",
                           values = APOE_Mouse$gene_id)

APOE_Mouse %<>%
  left_join(orthologs_mouse, by = c("gene_id" = "Gene.stable.ID")) %>%
  dplyr::select(-c(gene_id)) %>%
  dplyr::rename(gene_id  = HGNC.symbol) %>%
  group_by(gene_id) %>%
  dplyr::summarize_if(is.numeric, mean) %>%
  na.omit() %>%
  column_to_rownames(var = 'gene_id') %>%
  mutate_all( ~2^.x) %>% #undoinglogtransformation
  normalizeQuantiles() 

Mouse_obj <- CreateSeuratObject(counts = APOE_Mouse)

JakeLusis %<>%
  left_join(Strain_Batch, by = c("categories" = "Name"))

 
Mouse_obj@meta.data %<>%
  dplyr::mutate(
    Batch = Strain_Batch$Batch)  

Mouse_obj  %<>%
  FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData( verbose = FALSE, features = rownames(Mouse_obj)) %>%
  RunPCA( verbose = FALSE) %>%
  RunUMAP( dims = 1:12)

Umap_Before_Batch <- DimPlot(Mouse_obj, reduction = "umap",group.by = c('Batch')  , pt.size = 1.5) + labs(title = "Batch Before")
Umap_Before_Batch

PCA_Before_Batch <- DimPlot(Mouse_obj, reduction = "pca",group.by = c('Batch')  , pt.size = 1.5) + labs(title = "Batch Before")
PCA_Before_Batch

#Remove Batch  ----------------------------------------------------
Mouse_Obj.list <- SplitObject(Mouse_obj, split.by = "Batch")


Mouse_Obj.list <- lapply(X = Mouse_Obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nrow(x))
})

features_mouse <- SelectIntegrationFeatures(object.list = Mouse_Obj.list, nfeatures = 5000)

Mouse_Obj.anchors <- FindIntegrationAnchors(object.list = Mouse_Obj.list, anchor.features = features_mouse, k.filter = 50)
Mouse_Obj.combined <- IntegrateData(anchorset = Mouse_Obj.anchors, k.weight = 50)

DefaultAssay(Mouse_Obj.combined) <- "integrated"


saveRDS(Mouse_Obj.combined, file = "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/SeuratObjects/Jake_BatchCorrect.rds")


as.data.frame(GetAssayData(object = Mouse_Obj.combined, assay = "RNA", slot = "data")) -> Integrated_Mouse_all_genes
as.data.frame(GetAssayData(object = Mouse_Obj.combined, assay = "integrated", slot = "data")) -> Integrated_Mouse

#Integrated_Mouse
write.table((as.data.frame(GetAssayData(object = Mouse_Obj.combined, assay = "integrated", slot = "data"))%>% rownames_to_column(var = "gene_id")), 
            "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569_Corrected_Batch.csv", 
            sep = ',', row.names = F, col.names = T, quote = F)
#Integrated_Mouse_all_genes
write.table((as.data.frame(GetAssayData(object = Mouse_Obj.combined, assay = "RNA", slot = "data"))%>% rownames_to_column(var = "gene_id")), 
            "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569_Corrected_Batch_All_Genes.csv", 
            sep = ',', row.names = F, col.names = T, quote = F)

HumanPath <- "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Human/Human1/bulk_RNAseq_raw_counts.txt.minRib.txt.PC.txt"

Human_GeneSymb <- read_delim(HumanPath, show_col_types = FALSE) %>%
  dplyr::select(-c(start,end,ensembl_gene_id,  strand)) %>%
  dplyr::rename(gene_id = symbol) %>%
  remove_rownames() %>%
  group_by(gene_id) %>%
  dplyr::summarize_if(is.numeric, mean) %>%
  drop_na() %>%
  column_to_rownames(var = 'gene_id') %>%
  mutate_all( ~ifelse(.x > 4095, 4095, .x)) %>%
  mutate_all( ~-4096*log(1-(.x/4096))) %>%
  normalizeQuantiles()  %>%
  rownames_to_column(var = "gene_id") 

#write.table(Human_GeneSymb, file = paste0(file_path_sans_ext(HumanPath), "_", basename(dirname(HumanPath)),"Gene_Symbol" ,".txt"), sep = "\t", row.names = FALSE)
Human_GeneSymb <- read_delim(file = paste0(file_path_sans_ext(HumanPath), "_", basename(dirname(HumanPath)),"Gene_Symbol" ,".txt")) 

Integrated_Mouse %>%
  rownames_to_column(var = "gene_id") %>%
  left_join(Human_GeneSymb, by = c("gene_id")) %>%
  distinct(gene_id , .keep_all = TRUE) %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'gene_id') %>%
  normalizeQuantiles() -> Human_Mice

#for writing table!

Human_Mice %<>%
  rownames_to_column(var = 'gene_id')
#write.table(Human_Mice, file = paste0(file_path_sans_ext(HumanPath), "_", basename(dirname(HumanPath)),"Human_Mice" ,".txt"), sep = "\t", row.names = FALSE)
Human_Mice <- read_delim(file = paste0(file_path_sans_ext(HumanPath), "_", basename(dirname(HumanPath)),"Human_Mice" ,".txt")) 
