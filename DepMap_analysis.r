#Read file----
model <- read.csv2("~/Downloads/Model.csv", header = T, sep = ",")
gene_profile <- read.csv2("~/Downloads/OmicsExpressionGeneSetEnrichment.csv", header = T, sep = ",")
gsea <- read.csv2("~/Downloads/OmicsExpressionGeneSetEnrichment-2.csv", header = T, sep = ",")
somatic_mutations <- read.csv2("~/Downloads/OmicsSomaticMutations.csv", header = T, sep = ",")

#Data_Preparation----
model_sar <- model %>% dplyr::filter(OncotreeSubtype %in% types_of_sarcoma)
model_id_sar <- model_sar$ModelID
geneexp_sarc <- gene_profile %>% dplyr::filter(X %in% model_id_sar)
rownames(geneexp_sarc) <- geneexp_sarc$X
geneexp_sarc <- geneexp_sarc[,-1]
geneexp_sarc_I <- t(geneexp_sarc)
geneexp_sarc_I <- as.data.frame(geneexp_sarc_I)
geneexp_sarc_I <- as.data.frame(sapply(geneexp_sarc_I, as.numeric))
rownames(geneexp_sarc_I) <- colnames(geneexp_sarc)

rownames(model_sar)<-model_sar$ModelID
gene_model <- intersect(rownames(model_sar), colnames(geneexp_sarc_I))
model_sar <- model_sar[gene_model,]
geneexp_sarc_I <- geneexp_sarc_I[,gene_model]


samples <- colnames(geneexp_sarc_I)
condition_ts <- c()

for(i in 1:length(samples)){
  if(model_sar$OncotreeSubtype[i]=="Dedifferentiated Liposarcoma"){
    condition_ts <- c(condition_ts, "DDLPS")
  }
  if(model_sar$OncotreeSubtype[i]=="Leiomyosarcoma"){
    condition_ts <- c(condition_ts, "LMS")
  }
  if(model_sar$OncotreeSubtype[i]=="Synovial Sarcoma"){
    condition_ts <- c(condition_ts, "SS")
  }
  if(model_sar$OncotreeSubtype[i]=="Malignant Peripheral Nerve Sheath Tumor"){
    condition_ts <- c(condition_ts, "MPNST")
  }
  if(model_sar$OncotreeSubtype[i]=="Undifferentiated Pleomorphic Sarcoma/Malignant Fibrous Histiocytoma/High-Grade Spindle Cell Sarcoma"){
    condition_ts <- c(condition_ts, "UPS")
  }
}

model_sar$TS <- condition_ts

#Diner----

#PCA----
pca.res_dp <- pca(geneexp_sarc_I, metadata=model_sar)
biplot(pca.res_dp)
#pdf("PCA_dp.pdf")
biplot(pca.res_dp, colby = "TS", hline = 0, vline = 0,legendPosition = 'top',
       title = "PCA in cell lines") # Biplot with colors by sample type
#dev.off()
biplot(pca.res_dp, lab="", colby="OncotreeSubtype", hline = 0, vline = 0,legendPosition = 'top') # Biplot without sample names
biplot(pca.res_dp, x="PC1", y="PC3",lab="",colby="OncotreeSubtype", hline = 0, vline = 0,legendPosition = 'top') # Biplot with PC1 and PC3
pairsplot(pca.res_dp, colby="OncotreeSubtype")

#GSEA----
rownames(gsea)<-gsea$depmap_id

gsea_an <- gsea[,-c(1:8)]
gsea_an <- t(gsea_an)
gsea_an <- as.data.frame(gsea_an)

col_names_gsea_I <- sapply(colnames(gsea_an), function(x)unlist(strsplit(x,"\\-"))[1])
col_names_gsea_II <- sapply(colnames(gsea_an), function(x)unlist(strsplit(x,"\\-"))[2])
col_names_gsea <- paste(col_names_gsea_I, col_names_gsea_II, sep="")
colnames(gsea_an) <- col_names_gsea

write(cbind(gsea$depmap_id,gsea$lineage_3), "cell_id.txt")

#Genomics----
sm_sar <- somatic_mutations %>% filter(ModelID %in% model_id_sar)

View(sm_sar)