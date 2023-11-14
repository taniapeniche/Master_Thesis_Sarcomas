#Select the Human Data
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
#Create list of Kegg pathways
hs_kegg_df <- hs_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C2", # This is to filter only to the C2 curated gene sets
    gs_subcat == "CP:KEGG" # This is because we only want KEGG pathways
  )

#####Overall Analysis####
#Hub Genes in the blue module for the Leiomyosarcoma
#Gene List
blue_over <- blue$EST
#ORA Analysis
kegg_ora_results <- enricher(
  gene = blue_over , # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "fdr", # Method to be used for multiple testing correction
  universe = rownames(mRNA_fpkm), # A vector containing your background set genes
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_kegg_df,
    gs_name,
    gene_symbol
  )
)

#Signicant Results
kegg_result_df <- data.frame(kegg_ora_results@result)
sign_kegg_result <- kegg_result_df %>% filter(p.adjust <= 0.05)
#Plots
enrich_plot <- enrichplot::dotplot(kegg_ora_results)
cnetplot(kegg_ora_results, categorySize="pvalue")

#####Leyomiosarcomas####
#Hub Genes in the blue module for gender in leiomyosarcomas
#Gene List
blue_lms_genes<- blue_lms$EST
#ORA Analysis
kegg_ora_results_lms <- enricher(
  gene = blue_lms_genes , # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "fdr", # Method to be used for multiple testing correction
  universe = rownames(lms_mrna_fpkm), # A vector containing your background set genes
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_kegg_df,
    gs_name,
    gene_symbol
  )
)

#Signicant Results
kegg_result_df_lms <- data.frame(kegg_ora_results_lms@result)
sign_kegg_result_lms <- kegg_result_df_lms %>% filter(p.adjust <= 0.05)

#####MPNST####
#Hub Genes in MPNST
#Gene List
modules <- c("salmon4","black", "mb","dg","maroon")
mpnst_hgm <- vector("list",length = length(modules))
names(mpnst_hgm)<-modules

mpnst_hgm[[1]] <- salmon4_mpnst$EST
mpnst_hgm[[2]] <- black_mpnst$EST
mpnst_hgm[[3]] <- mb_mpnst$EST
mpnst_hgm[[4]] <- dg_mpnst$EST
mpnst_hgm[[5]] <- maroon_mpnst$EST

#ORA Analysis
kegg_ora_results_mpnst <- vector("list",length = length(modules))
names(kegg_ora_results_mpnst) <- modules

for(i in 1:length(kegg_ora_results_mpnst)){
  kos<- enricher(
    gene = mpnst_hgm[[i]] , # A vector of your genes of interest
    pvalueCutoff = 0.05, # Can choose a FDR cutoff
    pAdjustMethod = "fdr", # Method to be used for multiple testing correction
    universe = rownames(mpnst_mrna_fpkm), # A vector containing your background set genes
    # The pathway information should be a data frame with a term name or
    # identifier and the gene identifiers
    TERM2GENE = dplyr::select(
      hs_kegg_df,
      gs_name,
      gene_symbol
    )
  )
  kegg_ora_results_mpnst[[i]] <- kos
}

#Signicant Results
kegg_ora_df_mpnst <- vector("list",length = length(modules))
names(kegg_ora_df_mpnst) <- modules
sign_kegg_result_mpnst <- vector("list",length = length(modules))
names(sign_kegg_result_mpnst) <- modules

for(i in 1:length(kegg_ora_df_mpnst)){
  df <- data.frame(kegg_ora_results_mpnst[[1]]@result)
  kegg_ora_df_mpnst[[i]] <- df
  sr <- df %>% filter(p.adjust <= 0.05)
  sign_kegg_result_mpnst[[i]] <- sr
}

#####Synovial Sarcomas####
#Hub Genes in Synovial Sarcoma
#Gene List
modules_ss <- c("coral","ls", "dg")
ss_hgm <- vector("list",length = length(modules))
names(ss_hgm)<-modules

ss_hgm[[1]] <- coral_ss$EST
ss_hgm[[2]] <- ls_ss$EST
ss_hgm[[3]] <- dg_ss$EST

#ORA Analysis
kegg_ora_results_ss <- vector("list",length = length(modules))
names(kegg_ora_results_ss) <- modules

for(i in 1:length(kegg_ora_results_ss)){
  kos<- enricher(
    gene = ss_hgm[[i]] , # A vector of your genes of interest
    pvalueCutoff = 0.05, # Can choose a FDR cutoff
    pAdjustMethod = "fdr", # Method to be used for multiple testing correction
    universe = rownames(ss_mrna_fpkm), # A vector containing your background set genes
    # The pathway information should be a data frame with a term name or
    # identifier and the gene identifiers
    TERM2GENE = dplyr::select(
      hs_kegg_df,
      gs_name,
      gene_symbol
    )
  )
  kegg_ora_results_ss[[i]] <- kos
}

#Signicant Results
kegg_ora_df_ss <- vector("list",length = length(modules))
names(kegg_ora_df_ss) <- modules
sign_kegg_result_ss <- vector("list",length = length(modules))
names(sign_kegg_result_ss) <- modules

for(i in 1:length(kegg_ora_df_ss)){
  df <- data.frame(kegg_ora_results_ss[[1]]@result)
  kegg_ora_df_ss[[i]] <- df
  sr <- df %>% filter(p.adjust <= 0.05)
  sign_kegg_result_ss[[i]] <- sr
}




