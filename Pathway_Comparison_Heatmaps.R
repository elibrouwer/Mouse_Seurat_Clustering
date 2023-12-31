pacman::p_load("KEGGREST","stringr", "ComplexHeatmap", "circlize", "reshape", "reshape2")

#from https://gtexportal.org/home/datasets
Healthy_Aorta <- read.delim(file="C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/gene_tpm_2017-06-05_v8_artery_aorta.gct.gz", skip=2) %>%
  dplyr::select(-c(id,Name)) %>%
  dplyr::rename(gene_id = Description) %>%
  remove_rownames() %>%
  group_by(gene_id) %>%
  dplyr::summarize_if(is.numeric, mean) %>%
  column_to_rownames(var = 'gene_id') %>%
  drop_na() %>%
  normalizeQuantiles() %>%
  dplyr::select(sample(seq_len(ncol(.)), size = 150)) %>% #downsampling the aorta samples to n=150
  rownames_to_column(var = "gene_id") 


Human_Carotid_GeneSymb <- read_delim(file = paste0(file_path_sans_ext(HumanPath), "_", basename(dirname(HumanPath)),"Gene_Symbol" ,".txt")) #Human carotid samples with gene symbols

Integrated_Mouse_all_genes %>%
  rownames_to_column(var = "gene_id") %>%
  left_join(Human_Carotid_GeneSymb, by = c("gene_id")) %>%
  distinct(gene_id , .keep_all = TRUE) %>%
  left_join(Healthy_Aorta, by = c("gene_id")) %>%
  distinct(gene_id , .keep_all = TRUE) %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'gene_id') %>%
  normalizeQuantiles() -> Carotid_Mice_Aorta #was Human_Mice_Aorta

#write.table(Carotid_Mice_Aorta %>% rownames_to_column(var = 'gene_id'), file = paste0(file_path_sans_ext(MousePath), "_", basename(dirname(MousePath)),"Carotid_Mice_Aorta" ,".txt"), sep = "\t", row.names = FALSE)

Carotid_Mice_Aorta <- read_delim(file = paste0(file_path_sans_ext(MousePath), "_", basename(dirname(MousePath)),"Carotid_Mice_Aorta" ,".txt")) %>%
  column_to_rownames(var ='gene_id')

# Add data Coronary ----------------------------------------------------------

metadata_clint <- read_delim(file="C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Human/Human2/Metadata_Clint_Human2.txt") %>%
  dplyr::filter(Tissue =="Coronary_artery") %>%
  dplyr::select(c(SampleID, Classification)) 

metadata_clint %>%
  dplyr::filter(Classification == "Ischemic" | Classification == "Nonischemic") -> metadata_clint_Isch_vs_NonIsch


data_clint_no_aorta <-  read_delim(file="C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Human/Human2/genecountsraw_no_aorta.txt.PC.txt.minRib2.txt") %>%
  dplyr::select(-c(ENSID)) %>%
  dplyr::rename(gene_id = Description) %>%
  remove_rownames() %>%
  group_by(gene_id) %>%
  dplyr::summarize_if(is.numeric, mean) %>%
  drop_na() %>%
  column_to_rownames(var = 'gene_id') %>%
  normalizeQuantiles()  %>%
  rownames_to_column(var = "gene_id") 

#data_clint_no_aorta %<>%
#  column_to_rownames(var = "gene_id") %>%
#  dplyr::select(c(metadata_clint_Isch_vs_NonIsch$SampleID)) %>% #can be used to select only the ischemic vs non ischemic
#  rownames_to_column(var = "gene_id") 

Carotid_Mice_Aorta %>%
  rownames_to_column(var = 'gene_id') %>%
  left_join(data_clint_no_aorta, by = c("gene_id")) %>%
  distinct(gene_id , .keep_all = TRUE) %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'gene_id') %>%
  normalizeQuantiles() -> Carotid_Mice_Aorta_Coronary


# Pathway_Comparison_Heatmaps ---------------------------------------------

PathwayComparison_plt <- function(Mouse, Carotids, Coronary, Aorta, Selected_Paths, General_Paths){
  if(Mouse == TRUE & Carotids == TRUE & Coronary == FALSE & Aorta == FALSE){
    print("Mouse vs Caro")
    Integrated_Mouse_all_genes %>%
      rownames_to_column(var = "gene_id") %>%
      left_join(Human_GeneSymb, by = c("gene_id")) %>%
      distinct(gene_id , .keep_all = TRUE) %>%
      na.omit() %>%
      remove_rownames() %>%
      column_to_rownames(var = 'gene_id') %>%
      normalizeQuantiles() ->> Combined_Expression_dataframe
  } else if(Mouse == TRUE & Carotids == FALSE & Coronary == TRUE & Aorta == TRUE) {
    print("Mouse vs Cor vs Aorta")
    Integrated_Mouse_all_genes %>%
      rownames_to_column(var = "gene_id") %>%
      left_join(data_clint_no_aorta, by = c("gene_id")) %>%
      distinct(gene_id , .keep_all = TRUE) %>%
      left_join(Healthy_Aorta, by = c("gene_id")) %>%
      distinct(gene_id , .keep_all = TRUE) %>%
      na.omit() %>%
      remove_rownames() %>%
      column_to_rownames(var = 'gene_id') %>%
      normalizeQuantiles() ->> Combined_Expression_dataframe
    
  } else if(Mouse == FALSE & Carotids == TRUE & Coronary == TRUE & Aorta == TRUE) {
    print("HumanCompare")
    Human_GeneSymb %>%
      left_join(data_clint_no_aorta, by = c("gene_id")) %>%
      distinct(gene_id , .keep_all = TRUE) %>%
      left_join(Healthy_Aorta, by = c("gene_id")) %>%
      distinct(gene_id , .keep_all = TRUE) %>%
      na.omit() %>%
      remove_rownames() %>%
      column_to_rownames(var = 'gene_id') %>%
      normalizeQuantiles() ->> Combined_Expression_dataframe
  } else if(Mouse == TRUE & Carotids == TRUE & Coronary == TRUE & Aorta == TRUE) {
    print("all")
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
      normalizeQuantiles() ->> Combined_Expression_dataframe
  }
  Seurat_Obj <- CreateSeuratObject(counts = Combined_Expression_dataframe)
  
  Seurat_Obj@meta.data %<>%
    rownames_to_column("Names") %>%
    dplyr::mutate(dataset = ifelse(grepl("ae",Names), paste0("Human_Caro"), 
                                   ifelse(grepl("GSM",Names), paste0("Mouse"),
                                          ifelse(grepl("GTEX",Names), paste0("Human_Aorta"),
                                                 ifelse(grepl("UVA",Names), paste0("Human_Cor"), "-"))))) %>%
    column_to_rownames("Names")
  
  Seurat_Obj@meta.data %<>%
    rownames_to_column(var = 'Name') %>%
    left_join(Human_Mouse_Cats, by = c("Name" = "categories")) %>%
    left_join(SeuratClusters, by = c("Name" = "study_number")) %>%
    dplyr::mutate(cluster = coalesce(cluster, cats)) %>%
    mutate_at('cluster',  ~replace_na(., "Human_Aorta"))  %>%
    column_to_rownames("Name")
  
  Seurat_Obj %<>%
    FindVariableFeatures() %>%
    ScaleData(verbose = FALSE, features = rownames(Seurat_Obj)) 
  
  #calculate module scores
  for (j in 1: length(Selected_Paths$ID)){
    Path_Name <- query(Selected_Paths[j,])$displayName
    orthology <- getOrthology(Selected_Paths[j,], 'Homo sapiens')$stId
    print(Path_Name)
    print(Selected_Paths[j,])
    pathway=c(Selected_Paths[j,])
    path_genes <- getLDS(mart = Ens_Mouse,
                         attributes = 'reactome',
                         martL = Ens_Human,
                         attributesL = c("hgnc_symbol"),
                         filters = 'reactome',
                         values = pathway)
    print(path_genes$HGNC.symbol)
    Seurat_Obj %<>%
      AddModuleScore(features = list(c(path_genes$HGNC.symbol)) , name = paste0( make.names(Path_Name), "."), nbin = 24, ctrl = 100, assay = "RNA", slot = "scale.data", seed = 1)
    #print(VlnPlot(Human_Mouse_Obj.combined, features =paste0( make.names(Path_Name),".",  "1")))
    #print(DoHeatmap(Aorta_Human_Mouse_Obj, features = path_genes$HGNC.symbol, group.by = "cluster") + labs(title = paste0(Path_Name)))
  }

  #filtered based on signature that the pthws always ends with .1
  Names_Select <- colnames(Seurat_Obj@meta.data) %>%
    str_extract_all(pattern = ".+(?<!general)\\.1$") %>%
    unlist() %>%
    as.data.frame() 
  #modulescores are fetched from metadata and put in a new assay 'module'
  #data is centered and scaled
  Seurat_Obj[['module']] <- CreateAssayObject(data = t(x = FetchData(object = Seurat_Obj, vars = Names_Select$.)))
  Seurat_Obj %<>%
    ScaleData(assay = "module", do.scale = T, do.center = T)
  
 
  Humans_No_PlaqCat <- c("ae30", "ae245") #these humans dont have a plaque cat assigned so no cat, they are left out
  
  Seurat_Obj[['module']]@scale.data[,!colnames(Seurat_Obj) %in% Humans_No_PlaqCat ] %>%
    as.matrix()  -> mat
  
  rownames(mat) <- str_replace_all(str_sub_all(rownames(mat), end = -3 ) ,"\\."," ")
  
  mat <- mat[order(row.names(mat)),]
  
  cluster_anno <- Seurat_Obj[,!colnames(Seurat_Obj) %in% Humans_No_PlaqCat ]@meta.data$cluster
  rownames_anno <- rownames(mat)
  
  col_fun <- colorRamp2(c(-2,0,2), c("darkblue","grey90","darkorange"))
  
  selected_heatmap <<- Heatmap(mat, name = "Normalized Module Scores",  
                                               column_split = factor(cluster_anno),
                                               cluster_columns = F,
                                               show_column_dend = FALSE,
                                               cluster_column_slices = TRUE,
                                               column_title_gp = gpar(fontsize = 10),
                                               row_title_gp = gpar(fontsize = 10),
                                               column_gap = unit(1., "mm"),
                                               #column_title = "Atherosclerosis pathways which drive the clustering" ,
                                               cluster_rows = F,
                                               show_row_dend = FALSE,
                                               col = col_fun,
                                               row_names_gp = gpar(fontsize = 10),
                                               column_title_rot = 0,
                                               #top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
                                               left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#F8766D" , "#93AA00", "#00BA38", "#00C19F" ,"#00B9E3" ,"#619CFF" ,"#DB72FB" ,"#FF61C3")))),
                                               show_column_names = FALSE,
                                               use_raster = TRUE,
                                               raster_quality = 5,
                                               row_names_side = "right",
                                               row_title_side = "left",
                                               row_split = c("Mouse_2", "Mouse_1", "Mouse_1", "Mouse_2", "Mouse_4", 
                                                             "Mouse_3",  "Mouse_4", "Mouse_3"),
                                               row_title_rot = 0,
                                               column_names_rot = 45)
  
  #same but then for general paths
  for (j in 1: length(General_Paths$ID)){
    print(General_Paths[j,])
    pathway=c(General_Paths[j,])
    if (startsWith(pathway, "R")){
      Path_Name <- query(General_Paths[j,])$displayName 
      path_genes <- getLDS(mart = Ens_Mouse,
                           attributes = 'reactome',
                           martL = Ens_Human,
                           attributesL = c("hgnc_symbol"),
                           filters = 'reactome',
                           values = pathway)
      
      path_gen_list <- path_genes$HGNC.symbol
      print("reactome")
    } else {
      print("no reactome")
      Path_Name <- keggGet(pathway)[[1]]$NAME%>%
        str_replace("(- Mus musculus \\(house mouse\\))", "")
      keggGet(pathway)[[1]]$GENE -> names
      namesodd <-  names[seq(0,length(names),2)]
      path_gen_list <- gsub("\\;.*","",namesodd) %>%
        toupper()
    }
    print(Path_Name)
    print(path_gen_list)
    Seurat_Obj %<>%
      AddModuleScore(features = list(c(path_gen_list)) , name = paste0( make.names(Path_Name), ".general."), nbin = 24, ctrl = 100, assay = "RNA", slot = "scale.data", seed = 1)
    #print(VlnPlot(Human_MiceSeurat_Obj_Aorta_Coronary_Obj, features =paste0( make.names(Path_Name),".general.",  "1"), group.by = "cluster"))
    #print(DoHeatmap(Human_Mice_Aorta_Coronary_Obj, features = path_gen_list, group.by = "cluster", assay = "RNA", slot = "scale.data") + labs(title = paste0(Path_Name)))
  }
  
  Names_Select_Aorta_Coronary_General <- colnames(Seurat_Obj@meta.data) %>%
    str_extract_all(pattern = ".+.general\\.1$") %>%
    unlist() %>%
    as.data.frame() 
  
  Seurat_Obj[['module']] <- CreateAssayObject(data = t(x = FetchData(object = Seurat_Obj, vars = Names_Select_Aorta_Coronary_General$.)))
  Seurat_Obj %<>%
    ScaleData(assay = "module", do.scale = T, do.center = T, scale.max = 2)

  Humans_No_PlaqCat <- c("ae30", "ae245") #these humans dont have a plaque cat assigned so no cat, they are left out
  Seurat_Obj[['module']]@data -> a
  Seurat_Obj[['module']]@scale.data[,!colnames(Seurat_Obj) %in% Humans_No_PlaqCat ] %>%
    as.matrix()  -> mat_general
  rownames(mat_general) <- str_replace_all(str_sub_all(rownames(mat_general), end = -11 ) ,"\\."," ")
  cluster_anno_general <- Seurat_Obj[,!colnames(Seurat_Obj) %in% Humans_No_PlaqCat ]@meta.data$cluster

  
  General_paths_htmap <<- Heatmap(mat_general, name = "Normalized Module Scores",  
                                 column_split = factor(cluster_anno_general),
                                 cluster_columns = F,
                                 show_column_dend = FALSE,
                                 cluster_column_slices = TRUE,
                                 column_title_gp = gpar(fontsize = 10),
                                 row_title_gp = gpar(fontsize = 10),
                                 column_gap = unit(1., "mm"),
                                 cluster_rows = F,
                                 show_row_dend = FALSE,
                                 #column_title = "General atherosclerosis pathways", 
                                 col = col_fun,
                                 row_names_gp = gpar(fontsize = 10),
                                 column_title_rot = 0,
                                 top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
                                 left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#F8766D" , "#93AA00", "#00BA38", "#00C19F" ,"#00B9E3" ,"#619CFF" ,"#DB72FB" ,"#FF61C3")))),
                                 show_column_names = FALSE,
                                 use_raster = TRUE,
                                 row_split = c("Inflammatory","Metabolism" ,"Apoptosis", "ECM modification", "Angiogenesis", "Metabolism", "Angiogenesis", "Inflammatory", "Inflammatory", 
                                               "Inflammatory", "Inflammatory", "Shear stress"),
                                 raster_quality = 10,
                                 row_names_side = "right",
                                 row_title_side = "left",
                                 row_title_rot = 0)
  
}
#first made a function because the first idea was to make different graphs, now still using the function because I then save less big datasets in global environment.

PathwayComparison_plt(TRUE,TRUE,TRUE,TRUE, Selected_Paths, General_Paths ) #using the integrated coronary and carotid data.

General_Paths <- data.frame(ID = c( "R-MMU-6798695","R-MMU-70171", "R-MMU-2559580", "R-MMU-1474228", "R-MMU-194138", "mmu05417", "mmu04370", "mmu04064", "mmu04145", "mmu04668", "mmu04350", "mmu05418"))

Selected_Paths <- read.csv("C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569_Topn2paths.csv") 

paths <- c( "R-MMU-6798695", "R-MMU-71403", "R-MMU-70171", "R-MMU-445355" , "R-MMU-77289", "R-MMU-168256","R-MMU-6798695","R-MMU-2559580", "R-MMU-1474228", "R-MMU-194138", "mmu05417", "mmu04370", "mmu04064", "mmu04145", "mmu04668", "mmu04350", "mmu05418")

ht_list = General_paths_htmap %v% selected_heatmap #get the heatmaps vertical of eachother

draw(ht_list,column_title = "Atherosclerosis pathways", ht_gap = unit(.5, "cm"), heatmap_legend_side = "left",annotation_legend_side = "bottom", width = unit(40, "cm"))


# Pathway_Correlation_Heatmps ---------------------------------------------

Carotid_Mice_Aorta_Coronary_Obj <- CreateSeuratObject(counts = Combined_Expression_dataframe)
Carotid_Mice_Aorta_Coronary_Obj %<>%
  FindVariableFeatures() %>%
  ScaleData(verbose = FALSE, features = rownames(Seurat_Obj)) 

pdf(file= paste0(FolderPath, "Path_Correlation.pdf"), width = 4, height = 3)
for (j in 1:length(paths)){
  print(paths[j])
  pathway=c(paths[j])
if (startsWith(pathway, "R")){ #two different sources reactome (start with "R")
  Path_Name <- query(paths[j])$displayName 
  path_genes <- getLDS(mart = Ens_Mouse,
                       attributes = 'reactome',
                       martL = Ens_Human,
                       attributesL = c("hgnc_symbol"),
                       filters = 'reactome',
                       values = pathway)
  
  path_gen_list <- path_genes$HGNC.symbol
  print("reactome")
} else { #other database Kegg
  print("no reactome")
  Path_Name <- keggGet(pathway)[[1]]$NAME%>%
    str_replace("(- Mus musculus \\(house mouse\\))", "")
  keggGet(pathway)[[1]]$GENE -> names
  namesodd <-  names[seq(0,length(names),2)]
  path_gen_list <- gsub("\\;.*","",namesodd) %>% #have to select all the odd names because output of keggGet is "gene" and number
    toupper()
}

GetAssayData(Carotid_Mice_Aorta_Coronary_Obj, slot = "data") -> data
data %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  subset(gene_id %in% path_gen_list) %>%
  remove_rownames() %>%
  column_to_rownames(var = "gene_id") %>%
  dplyr::select(starts_with("ae")| starts_with("GSM")) %>% #select carotid (start with "ae") and Mouse (start with "GSM")
  cor() %>% #pearson's correlation
  data.frame() %>%
  rownames_to_column(var = "Name") %>%
  left_join(Human_Mouse_Cats, by =c("Name" = "categories")) %>%
  column_to_rownames(var = "Name") %>%
  group_by(cluster) %>%
  reshape2::melt() %>%
  left_join(Human_Mouse_Cats, by =c("variable" = "categories") ) %>%
  dplyr::filter(!grepl('Carotid_NA', c(cluster.x))) %>%
  dplyr::filter(!grepl('Carotid_NA', c(cluster.y))) %>% #leaving out the human_NA's 
  dplyr::filter(cluster.x >= cluster.y) %>% #filter to only be left with the lower triangle correlation matrix
   ggplot(aes(x = cluster.x, y = cluster.y, fill = value)) + geom_tile(color = "black", size=0.2) + ggtitle(paste0(Path_Name)) +   scale_fill_gradient2(low = "#075AFF",
                                                                                                                              mid = "white",
                                                                                                                              high = "darkgreen") + 
  theme_minimal()+
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_blank())+
  theme(axis.text.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 50), angle =0),
        axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 50), angle =90),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) -> correlation_plt

print(correlation_plt)
}

dev.list()
dev.off()


