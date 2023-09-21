#controlling the distribution of all the datasets
Human_Cor %>%
  column_to_rownames(var = 'gene_id') %>%
  melt() %>%
  dplyr::mutate(dataset = paste0("Human_Cor")) -> Human_Cor_mlt

Human_Mice_ldlr %>%
  melt() %>%
  dplyr::mutate(dataset = paste0("Mouse_ldl")) -> Mouse_ldl_mlt

Healthy_Aorta %>%
  column_to_rownames(var = 'gene_id') %>%
  melt() %>%
  dplyr::mutate(dataset = paste0("Human_Aorta")) -> Healthy_Aorta_mlt

Integrated_Mouse_all_genes %>%
  melt() %>%
  dplyr::mutate(dataset = paste0("Mouse")) -> Mouse_mlt

#control steps for the distribution
Human_Mice_Aorta_Coronary %>%
  melt() %>%
  dplyr::mutate(dataset = ifelse(grepl("ae",variable), paste0("Human_Caro"), 
                                 ifelse(grepl("GSM",variable), paste0("Mouse"),
                                        ifelse(grepl("GTEX",variable), paste0("Human_Aorta"),
                                               ifelse(grepl("UVA",variable), paste0("Human_Cor"), "-"))))) %>%
  ggplot(aes( x= value, color = dataset)) +
  geom_density(alpha=0.3,size=1) + scale_x_log10() + geom_density() + ggtitle("Density plot of all dataset after quantile-normalization") -> density_plt_after_normalization

Human_Caro %>%
  column_to_rownames(var = 'gene_id') %>%
  melt() %>%
  dplyr::mutate(dataset = paste0("Human_Caro")) %>%
  rbind(Human_Cor_mlt) %>%
  rbind(Mouse_mlt) %>%
  rbind(Healthy_Aorta_mlt) %>%
  #rbind(Mouse_ldl_mlt) %>%
  ggplot(aes( x= value, color = dataset)) +
  geom_density(alpha=0.3,size=1) + scale_x_log10() + geom_density() + ggtitle("Density plot of all dataset before quantile-normalization") -> density_plt_before_normalization