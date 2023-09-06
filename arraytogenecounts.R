library(biomaRt)
library(dplyr)
library(affy)
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
library(affyio)
library(oligo)
library(stringr)
library(data.table)
celpath = "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569" 
celFiles <- affy::list.celfiles(celpath)
data = ReadAffy(celfile.path = celpath) ## all the affyfiles are opened in celpath
cel <- read.celfiles(paste0(celpath, "/", celFiles))
data = oligo::rma(cel)

id <- featureNames(data)
annotationlist <- AnnotationDbi::select(htmg430pm.db, c(id), c("ENSEMBL"))

library(htmg430pm.db)


expressiondata <- data.frame(exprs(data)) ## expressiondataset from data

expressiondata %>%
  rownames_to_column(var = "PROBEID") %>%
  left_join(annotationlist, by = c("PROBEID")) %>%
  dplyr::select(-c(PROBEID)) %>%
  dplyr::rename(gene_id = ENSEMBL) %>%
  na.omit() %>%
  group_by(gene_id) %>%
  dplyr::summarize_if(is.numeric, mean)  -> annotated_expr

write.table(annotated_expr, "C:/Users/brouw/OneDrive - Universiteit Utrecht/Master Bioinformatics/Major Internship/R/Datasets/Mice/GSE66569_APOE/GSE66569.txt", sep = "\t", row.names= FALSE)

