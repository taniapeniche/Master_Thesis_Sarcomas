#TCGA----
#Summary_Plot-----
plotmafSummary(tcga_all_maf)

#Oncoplot-----
#genes
oncoplot(maf = tcga_all_maf, top = 10, sortByMutation = TRUE)
#Pathway
oncoplot(maf = tcga_all_maf, pathways = "auto", gene_mar = 8, fontSize = 0.6)

#Oncogenic Pathways-----
OncogenicPathways(maf = tcga_all_maf)

#Somatic Interactions-----
somaticInteractions(maf = tcga_all_maf, pvalue = c(0.05, 0.1), top = 20)

#Mutational Signatures-----
TCGA.tnm <- trinucleotideMatrix(maf = tcga_all_maf, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

##Signature analysis
TCGA.sign <- estimateSignatures(mat = TCGA.tnm, nTry = 6)

TCGA.sig <- extractSignatures(mat = TCGA.tnm, n = 5)

#Compate against original 30 signatures 
TCGA.og30.cosm <- compareSignatures(nmfRes = TCGA.sig, sig_db = "legacy")

#Compate against updated version3 60 signatures 
TCGA.v3.cosm <- compareSignatures(nmfRes = TCGA.sig, sig_db = "SBS")

#comparison of similarities of detected signatures against validated signatures
pheatmap(mat = TCGA.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

#plot signatures
plotSignatures(nmfRes = TCGA.sig, title_size = 1.2, sig_db = "SBS")

#Clinical Enrichment-----
TCGA.ce <- clinicalEnrichment(maf = tcga_all_maf, clinicalFeature = 'TUMOR_TISSUE_SITE')

#Results are returned as a list. Significant associations p-value < 0.05
TCGA.ce$groupwise_comparision[p_value < 0.05]

plotEnrichmentResults(enrich_res = TCGA.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)

#MSKCC------
#maf_summary_plot-----
plotmafSummary(maf = mskcc_all_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#Oncoplot------
oncoplot(maf = mskcc_all_maf, top = 10)
#Pathways
oncoplot(maf = mskcc_all_maf, pathways = "auto", gene_mar = 8, fontSize = 0.6)

#Somatic Interactions----
somaticInteractions(maf = MSK, top = 25, pvalue = c(0.05, 0.1))

#Mutational Signatures-----
MSK.tnm <- trinucleotideMatrix(maf = mskcc_all_maf, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

##Signature analysis
MSK.sign <- estimateSignatures(mat = MSK.tnm, nTry = 6)

MSK.sig <- extractSignatures(mat = MSK.tnm, n = 2)

#Compate against original 30 signatures 
MSK.og30.cosm <- compareSignatures(nmfRes = MSK.sig, sig_db = "legacy")

#Compate against updated version3 60 signatures 
MSK.v3.cosm <- compareSignatures(nmfRes = MSK.sig, sig_db = "SBS")

#comparison of similarities of detected signatures against validated signatures
pheatmap(mat = MSK.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

#plot signatures
plotSignatures(nmfRes = MSK.sig, title_size = 1.2, sig_db = "SBS")

#Clinical Enrichment------
#Primary Site
MSK_s.ce <- clinicalEnrichment(maf = mskcc_all_maf, clinicalFeature = 'PRIMARY_SITE')

#Results are returned as a list. Significant associations p-value < 0.05
MSK_s.ce$groupwise_comparision[p_value < 0.05]

plotEnrichmentResults(enrich_res = MSK_s.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)

#ASP_18-------
#maf_summary_plot-----
plotmafSummary(maf = asp_18, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#oncoplot-----
oncoplot(maf = asp_18, top = 20, sortByMutation = TRUE)
#Pathways
oncoplot(maf = asp_18, pathways = "auto", gene_mar = 8, fontSize = 0.6)

#Somatic Interactions----
somaticInteractions(maf = asp_18, top = 25, pvalue = c(0.05, 0.1))

#Compate against original 30 signatures
asp_18.og30.cosm <- compareSignatures(nmfRes = asp_18.sig, sig_db = "legacy")

#Compate against updated version3 60 signatures------
asp_18.v3.cosm <- compareSignatures(nmfRes = asp_18.sig, sig_db = "SBS")

#comparison of similarities of detected signatures against validated signatures
pheatmap(mat = asp_18.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

#plot signatures------
plotSignatures(nmfRes = asp_18.sig, title_size = 1.2, sig_db = "SBS")

#ASP_20-------
#maf_summary_plot-------
plotmafSummary(maf = maf_20, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

mafbarplot(maf = maf_20)

#oncoplot-------
oncoplot(maf = maf_20, top = 20, sortByMutation = TRUE)
#Pathways
oncoplot(maf = maf_20, pathways = "auto", gene_mar = 8, fontSize = 0.6)

#somatic Interactions----
somaticInteractions(maf = maf_20, pvalue = c(0.05, 0.1), top = 50)

#Mutational Signatures----
laml.tnm <- trinucleotideMatrix(maf = maf_20, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

##Signature analysis
laml.sign <- estimateSignatures(mat = laml.tnm, nTry = 12)

laml.sig <- extractSignatures(mat = laml.tnm, n = 8)

#Compate against original 30 signatures 
laml.og30.cosm <- compareSignatures(nmfRes = laml.sig, sig_db = "legacy")

#Compate against updated version3 60 signatures 
laml.v3.cosm <- compareSignatures(nmfRes = laml.sig, sig_db = "SBS")

#comparison of similarities of detected signatures against validated signatures
pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

#plot signatures
plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "SBS")

#####Associated with clinical feature####
fab.ce <- clinicalEnrichment(maf = maf_20, clinicalFeature = 'BX_LOCATION')

#Results are returned as a list. Significant associations p-value < 0.05
fab.ce$groupwise_comparision[p_value < 0.05]

plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)

#MAF Analyse per type of sarcoma----
#Select the samples of inteterest
all_maf_in<-all_maf[-c(5,10:13,16:19,22,23,25:32,35,37:40,42,45,47,53:56)]
#Save the name of the ones that  will be used
write(names(all_maf_in), "Types_of_Sarcoma_in_I.txt")
#DGEA
for(i in 1:length(all_maf_in)){
  maf <- all_maf_in[[i]]
  plotmafSummary(maf)
  #####Oncoplot####
  #genes
  oncoplot(maf = maf, top = 10, sortByMutation = TRUE, fontSize = 0.7)
  #Pathway
  oncoplot(maf = maf, pathways = "auto", gene_mar = 8, fontSize = 0.7)
  #####Oncogenic Pathways####
  OncogenicPathways(maf = maf)
}
#Somatic interections
for(i in 1:length(all_maf_in)){
  maf <- all_maf_in[[i]]
  #####Somatic Interactions####
  somaticInteractions(maf = maf, pvalue = c(0.05, 0.1), top = 20, fontSize = 0.7)
}

