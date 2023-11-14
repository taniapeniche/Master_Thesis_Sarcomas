#Leiomiosarcomas----
#Subbset for only 
both <- intersect(rownames(CI_TCGA_mRNA),names(ts_geneexp[["LMS"]]))

lms_mrna_fpkm <- mRNA_fpkm %>% 
  dplyr::select(all_of(both))

#Ci for ts
ci_tcga_lms <- CI_TCGA_mRNA %>% filter(Tumor_Sample_Barcode %in% both)

# Prepare gene expression matrix for PCA analysis
logfpkm_lms <- log2(lms_mrna_fpkm+1)
row.names(logfpkm_lms) <- row.names(lms_mrna_fpkm)

#Outliers 
gsg_lms <- goodSamplesGenes(t(logfpkm_lms))
summary(gsg_lms)
gsg_lms$allOK

table(gsg_lms$goodGenes)
table(gsg_lms$goodSamples)

# remove genes that are detectd as outliers
data_lms <- logfpkm_lms[gsg_lms$goodGenes == TRUE,]


# detect outlier samples - hierarchical clustering - method 1
htree_lms <- hclust(dist(t(data_lms)), method = "ward.D2")
plot(htree_lms)

#Remove outliers
data_lms <- data_lms %>% dplyr::select(!(c("TCGA-DX-A8BZ", "TCGA-JV-A75J", "TCGA-HS-A5N9")))
ci_tes<-ci_tcga_lms %>% filter(!(Tumor_Sample_Barcode %in% c("TCGA-DX-A8BZ", "TCGA-JV-A75J", "TCGA-HS-A5N9")))

#PCA
pca.res_lms <- pca(data_lms, metadata=ci_tes)

#Plot variance explained by each component
screeplot(pca.res_lms)

#Plot 2 selected components/eigenvectors
biplot(pca.res_lms)
biplot(pca.res_lms, hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
biplot(pca.res_lms, lab="", hline = 0, vline = 0,legendPosition = 'top') # Biplot without sample names
biplot(pca.res_lms, x="PC1", y="PC3",lab="",colby="TS", hline = 0, vline = 0,legendPosition = 'top') # Biplot with PC1 and PC3


#Plot several components
pairsplot(pca.res_lms)


#prepare date for dge
round_data_lms <- round(data_lms)

col_metada_lms <- ci_tes[,c(7,8,10:18,32,34,46:51,67:69)]

# create dds
dds_lms <- DESeqDataSetFromMatrix(countData = round_data_lms,
                                  colData = col_metada_lms,
                                  design = ~ 1)

## remove all genes with low counts
## suggested by WGCNA on RNAseq FAQ

dds75_lms <- dds_lms[rowSums(counts(dds_lms) >= 3) >= 24,]
nrow(dds75_lms)


# perform variance stabilization
dds_norm_lms <- varianceStabilizingTransformation(dds75_lms)

# get normalized counts
norm.counts_lms <- assay(dds_norm_lms) %>% t()


# Network Construction 
# Choose a set of soft-thresholding powers
power_lms <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function

sft_lms <- pickSoftThreshold(norm.counts_lms,
                             powerVector = power_lms,
                             networkType = "signed",
                             verbose = 5)

sft.data_lms <- sft_lms$fitIndices

# visualization to pick power

a1_lms <- ggplot(sft.data_lms, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2_lms <- ggplot(sft.data_lms, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1_lms, a2_lms, nrow = 2)


# convert matrix to numeric
norm.counts_lms[] <- sapply(norm.counts_lms, as.numeric)

soft_power_lms <- 8
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet_lms <- blockwiseModules(norm.counts_lms,
                              maxBlockSize = 14000,
                              TOMType = "signed",
                              power = soft_power_lms,
                              mergeCutHeight = 0.25,
                              numericLabels = FALSE,
                              randomSeed = 1234,
                              verbose = 3)


cor <- temp_cor

# Module Eigengenes 
module_eigengenes_lms <- bwnet_lms$MEs


# Print out a preview
head(module_eigengenes_lms)


# get number of genes for each module
table(bwnet_lms$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_lms$dendrograms[[1]], cbind(bwnet_lms$unmergedColors, bwnet_lms$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

moduleColors_lms <- bwnet_lms$colors


# grey module = all genes that doesn't fall into other modules were assigned to the grey module

# Relate modules to clinical traits
# module trait associations



# create traits file - binarize categorical variables
traits_lms <- col_metada_lms %>% 
  mutate(metastatisi = ifelse(grepl('YES', metastatic_disease_confirmed), 1, 0)) %>%
  mutate(gender = ifelse(grepl('FEMALE', gender), 1, 0)) %>%
  dplyr::select(1,3,6,23)

traits_lms <- traits_lms %>% mutate(status = ifelse(grepl('Alive', vital_status), 1, 0)) %>%
  mutate(otm = ifelse(grepl('Dead', history_other_malignancy), 1, 0))

traits_lms <- traits_lms %>% dplyr::select(1,4:6)


# binarize categorical variables


metastasis_site <- binarizeCategoricalColumns(col_metada_lms$metastasis_site,
                                              includePairwise = FALSE,
                                              includeLevelVsAll = TRUE,
                                              minCount = 1)

ethnicity <- binarizeCategoricalColumns(col_metada_lms$ethnicity,
                                   includePairwise = FALSE,
                                   includeLevelVsAll = TRUE,
                                   minCount = 1)

traits_lms <- cbind(traits_lms, metastasis_site,ethnicity)



# Define numbers of genes and samples
nSamples_lms <- nrow(norm.counts_lms)
nGenes_lms <- ncol(norm.counts_lms)


module.trait.corr_lms <- cor(module_eigengenes_lms, traits_lms, use = 'p')
module.trait.corr.pvals_lms <- corPvalueStudent(module.trait.corr_lms, nSamples_lms)



# visualize module-trait association as a heatmap

heatmap.data_lms <- merge(module_eigengenes_lms, traits_lms, by = 'row.names')

head(heatmap.data_lms)

heatmap.data_lms <- heatmap.data_lms[,-1]
colnames(heatmap.data_lms)[12:14] <- c( "Gender","Metastasis","Status")
colnames(heatmap.data_lms)[17:20] <- c("Bone", "Liver", "Lung", "Other") 
colnames(heatmap.data_lms)[22:24] <- c("Asian", "Black", "White")

pdf("lms_wgcna.pdf")
CorLevelPlot(heatmap.data_lms,
             x = names(heatmap.data_lms)[c(12:14)],
             y = names(heatmap.data_lms)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"),
             main = "WGCNA for clinical traits in LMS",
             cexMain = 1)
dev.off()


module.gene.mapping <- as.data.frame(bwnet_lms$colors)
module.gene.mapping %>% 
  filter(`bwnet_lms$colors` == 'blue') %>% 
  rownames()


#####Indentifying gene Hubs#
lms_g <- as.data.frame(traits_lms$gender)
names(lms_g) <- "lms_g"

#names (colors) of the modules
modNames_lms <-  substring(names(module_eigengenes_lms), 3)

geneModuleMembership_lms <- as.data.frame(cor(norm.counts_lms, module_eigengenes_lms, use = "p"))
MMPvalue_lms <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_lms), nSamples_lms))


names(geneModuleMembership_lms) <- paste("MM", modNames_lms, sep="")
names(MMPvalue_lms) <- paste("p.MM", modNames_lms, sep="")

geneTraitSignificance_lms <- as.data.frame(cor(norm.counts_lms, lms_g, use = "p"))
GSPvalue_lms <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_lms), nSamples_lms))
GSPvalue_lms_adjust <- as.data.frame(p.adjust(as.matrix(GSPvalue_lms), method = "fdr"))

names(geneTraitSignificance_lms) <- paste("GS.", names(lms_g), sep="")
names(GSPvalue_lms) <- paste("p.GS.", names(lms_g), sep="")
names(GSPvalue_lms_adjust) <- paste("p.adjust.GS.", names(lms_g), sep="")

module_lms <- "blue" #########################putting the color below the plot
column_lms <- match(module_lms, modNames_lms)
moduleGenes_lms <- moduleColors_lms==module_lms

verboseScatterplot(abs(geneModuleMembership_lms[moduleGenes_lms, column_lms]),
                   abs(geneTraitSignificance_lms[moduleGenes_lms, 1]),
                   xlab = paste("Module Membership in", module_lms, "module"),
                   ylab = "Gene significance for LMS",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_lms)

geneInfo0_lms <- data.frame(EST = colnames(norm.counts_lms),
                           moduleColor = moduleColors_lms,
                           geneTraitSignificance_lms,
                           GSPvalue_lms, GSPvalue_lms_adjust)


modOrder_lms <- order(-abs(cor(module_eigengenes_lms, lms_g, use = "p")))
for (mod in 1:ncol(geneModuleMembership_lms))
{
  oldNames = names(geneInfo0_lms)
  geneInfo0_lms = data.frame(geneInfo0_lms, geneModuleMembership_lms[, modOrder_lms[mod]], 
                             MMPvalue_lms[, modOrder_lms[mod]]);
  names(geneInfo0_lms) = c(oldNames, paste("MM.", modNames_lms[modOrder_lms[mod]], sep=""),
                           paste("p.MM.", modNames_lms[modOrder_lms[mod]], sep=""))
}
geneOrder_lms <- order(geneInfo0_lms$moduleColor, -abs(geneInfo0_lms$GS.lms_g))
geneInfo_lms <- geneInfo0_lms[geneOrder_lms, ]

blue_lms <- geneInfo_lms %>%
  filter(moduleColor=="blue") %>%
  filter(abs(GS.lms_g) >= 0.2) %>%
  filter(p.adjust.GS.lms_g <= 0.05) %>%
  filter(abs(MM.blue) >= 0.8) %>%
  filter(p.MM.blue <= 0.05)


#####Liposarcoma####
#Subbset for only ts
both <- intersect(rownames(CI_TCGA_mRNA),names(ts_geneexp[["DL"]]))

dl_mrna_fpkm <- mRNA_fpkm %>% 
  dplyr::select(all_of(both))

#Ci for ts
ci_tcga_dl <- CI_TCGA_mRNA %>% filter(Tumor_Sample_Barcode %in% both)

# Prepare gene expression matrix for PCA analysis: Step 1) get log TPMs (we add 1 TPM to each gene to avoid infinite values after log)
logfpkm_dl <- log2(dl_mrna_fpkm+1)
row.names(logfpkm_dl) <- row.names(dl_mrna_fpkm)

#Outliers 
gsg_dl <- goodSamplesGenes(t(logfpkm_dl))
summary(gsg_dl)
gsg_dl$allOK

table(gsg_dl$goodGenes)
table(gsg_dl$goodSamples)

# remove genes that are detectd as outliers
data_dl <- logfpkm_dl[gsg_dl$goodGenes == TRUE,]


# detect outlier samples - hierarchical clustering - method 1
htree_dl <- hclust(dist(t(data_dl)), method = "ward.D2")
plot(htree_dl)


#data_dl <- data_dl %>% select(!(c("TCGA-DX-A2IZ","TCGA-3B-A9HJ","TCGA-3R-A8YX")))
#ci_tcga_dl <- ci_tcga_dl %>% filter(!(Tumor_Sample_Barcode %in% c("TCGA-DX-A2IZ","TCGA-3B-A9HJ","TCGA-3R-A8YX")))
#PCA
pca.res_dl <- pca(data_dl, metadata=ci_tcga_dl)

#Plot variance explained by each component
screeplot(pca.res_dl)

#Plot 2 selected components/eigenvectors
biplot(pca.res_dl)
biplot(pca.res_dl, hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
biplot(pca.res_dl, lab="", hline = 0, vline = 0,legendPosition = 'top') # Biplot without sample names
biplot(pca.res_dl, x="PC1", y="PC3",lab="",colby="TS", hline = 0, vline = 0,legendPosition = 'top') # Biplot with PC1 and PC3


#Plot several components
pairsplot(pca.res_dl)

#Prepare Data for next step
round_data_dl <- round(data_dl)

col_metada_dl <- ci_tcga_dl[,c(7,8,10:18,32,34,46:51,67:69)]

# create dds
dds_dl <- DESeqDataSetFromMatrix(countData = round_data_dl,
                                  colData = col_metada_dl,
                                  design = ~ 1)

## remove all genes with low counts in all samples
## suggested by WGCNA on RNAseq FAQ

dds75_dl <- dds_dl[rowSums(counts(dds_dl) >= 3) >= 24,]
nrow(dds75_dl)


# perform variance stabilization
dds_norm_dl <- varianceStabilizingTransformation(dds75_dl)

# get normalized counts
norm.counts_dl <- assay(dds_norm_dl) %>% t()


# Network Construction 
# Choose a set of soft-thresholding powers
power_dl <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function

sft_dl <- pickSoftThreshold(norm.counts_dl,
                             powerVector = power_dl,
                             networkType = "signed",
                             verbose = 5)

sft.data_dl <- sft_dl$fitIndices

# visualization to pick power

a1_dl <- ggplot(sft.data_dl, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2_dl <- ggplot(sft.data_dl, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1_dl, a2_dl, nrow = 2)


# convert matrix to numeric
norm.counts_dl[] <- sapply(norm.counts_dl, as.numeric)

soft_power_dl <- 8
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet_dl <- blockwiseModules(norm.counts_dl,
                              maxBlockSize = 14000,
                              TOMType = "signed",
                              power = soft_power_dl,
                              mergeCutHeight = 0.25,
                              numericLabels = FALSE,
                              randomSeed = 1234,
                              verbose = 3)


cor <- temp_cor

# 5. Module Eigengenes 
module_eigengenes_dl <- bwnet_dl$MEs


# Print out a preview
head(module_eigengenes_dl)


# get number of genes for each module
table(bwnet_dl$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_dl$dendrograms[[1]], cbind(bwnet_dl$unmergedColors, bwnet_dl$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


moduleColors_dl <- bwnet_dl$colors

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

# Relate modules to clinical traits
# module trait associations



# create traits file - binarize categorical variables
traits_dl <- col_metada_dl %>% 
  mutate(metastatisi = ifelse(grepl('YES', metastatic_disease_confirmed), 1, 0)) %>%
  mutate(gender = ifelse(grepl('FEMALE', gender), 1, 0)) %>%
  dplyr::select(1,3,6,23)

traits_dl <- traits_dl %>% mutate(status = ifelse(grepl('Alive', vital_status), 1, 0)) %>%
  mutate(otm = ifelse(grepl('Dead', history_other_malignancy), 1, 0))

traits_dl <- traits_dl %>% dplyr::select(1,4:6)


# binarize categorical variables


metastasis_site <- binarizeCategoricalColumns(col_metada_dl$metastasis_site,
                                              includePairwise = FALSE,
                                              includeLevelVsAll = TRUE,
                                              minCount = 1)

ethnicity <- binarizeCategoricalColumns(col_metada_dl$ethnicity,
                                   includePairwise = FALSE,
                                   includeLevelVsAll = TRUE,
                                   minCount = 1)

traits_dl <- cbind(traits_dl, metastasis_site,ethnicity)


# Define numbers of genes and samples
nSamples_dl <- nrow(norm.counts_dl)
nGenes_dl <- ncol(norm.counts_dl)


module.trait.corr_dl <- cor(module_eigengenes_dl, traits_dl, use = 'p')
module.trait.corr.pvals_dl <- corPvalueStudent(module.trait.corr_dl, nSamples_dl)



# visualize module-trait association as a heatmap

heatmap.data_dl <- merge(module_eigengenes_dl, traits_dl, by = 'row.names')

head(heatmap.data_dl)

heatmap.data_dl <- heatmap.data_dl[,-1]
colnames(heatmap.data_dl)[9:11] <- c("Gender","Metastasis","Status")
colnames(heatmap.data_dl)[13:15] <- c("Lung", "Asian","White") 


pdf("dl_wgcna.pdf")
CorLevelPlot(heatmap.data_dl,
             x = names(heatmap.data_dl)[c(9:11,13:15)],
             y = names(heatmap.data_dl)[1:8],
             col = c("blue1", "skyblue", "white", "pink", "red"))

dev.off()

#see genes of the module
module.gene.mapping <- as.data.frame(bwnet_dl$colors)
module.gene.mapping %>% 
  filter(`bwnet_dl$colors` == 'tan') %>% 
  rownames()

#Indetify hub genes 
dl_a <- as.data.frame(traits_dl$Asian)
names(dl_a) <- "dl_a"

#names (colors) of the modules
modNames_dl <-  substring(names(module_eigengenes_dl), 3)

geneModuleMembership_dl <- as.data.frame(cor(norm.counts_dl, module_eigengenes_dl, use = "p"))
MMPvalue_dl <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_dl), nSamples_dl))


names(geneModuleMembership_dl) <- paste("MM", modNames_dl, sep="")
names(MMPvalue_dl) <- paste("p.MM", modNames_dl, sep="")

geneTraitSignificance_dl <- as.data.frame(cor(norm.counts_dl, dl_a, use = "p"))
GSPvalue_dl <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_dl), nSamples_dl))
GSPvalue_dl_adjust <- as.data.frame(p.adjust(as.matrix(GSPvalue_dl), method = "fdr"))

names(geneTraitSignificance_dl) <- paste("GS.", names(dl_a), sep="")
names(GSPvalue_dl) <- paste("p.GS.", names(dl_a), sep="")
names(GSPvalue_dl_adjust) <- paste("p.Adjust.GS.", names(dl_a), sep="")

module_dl <- "red" #########################putting the color below the plot
column_dl <- match(module_dl, modNames_dl)
moduleGenes_dl <- moduleColors_dl==module_dl

verboseScatterplot(abs(geneModuleMembership_dl[moduleGenes_dl, column_dl]),
                   abs(geneTraitSignificance_dl[moduleGenes_dl, 1]),
                   xlab = paste("Module Membership in", module_dl, "module"),
                   ylab = "Gene significance for dl",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_dl)

#Table with all the info
geneInfo0_dl <- data.frame(EST = colnames(norm.counts_dl),
                          moduleColor = moduleColors_dl,
                          geneTraitSignificance_dl,
                          GSPvalue_dl, GSPvalue_dl_adjust)


modOrder_dl <- order(-abs(cor(module_eigengenes_dl, dl_a, use = "p")))
for (mod in 1:ncol(geneModuleMembership_dl))
{
  oldNames = names(geneInfo0_dl)
  geneInfo0_dl = data.frame(geneInfo0_dl, geneModuleMembership_dl[, modOrder_dl[mod]], 
                            MMPvalue_dl[, modOrder_dl[mod]]);
  names(geneInfo0_dl) = c(oldNames, paste("MM.", modNames_dl[modOrder_dl[mod]], sep=""),
                          paste("p.MM.", modNames_dl[modOrder_dl[mod]], sep=""))
}
geneOrder_dl <- order(geneInfo0_dl$moduleColor, -abs(geneInfo0_dl$GS.dl_a))
geneInfo_dl <- geneInfo0_dl[geneOrder_dl, ]

#Selectc hub genes
red_dl <- geneInfo_dl %>%
  filter(moduleColor=="red") %>%
  filter(abs(GS.dl_a) >= 0.2) %>%
  filter(p.Adjust.GS.dl_a <= 0.05) %>%
  filter(abs(MM.green) >= 0.8) %>%
  filter(p.MM.green<= 0.05)

green_dl <- geneInfo_dl %>%
  filter(moduleColor=="green") %>%
  filter(abs(GS.dl_a) >= 0.2) %>%
  filter(p.Adjust.GS.dl_a <= 0.05) %>%
  filter(abs(MM.green) >= 0.8) %>%
  filter(p.MM.green<= 0.05)


#####Myxofibrossarcoma####
#Subbset for only ts
both <- intersect(rownames(CI_TCGA_mRNA),names(ts_geneexp[["MF"]]))

mf_mrna_fpkm <- mRNA_fpkm %>% 
  dplyr::select(all_of(both))

#Ci for ts
ci_tcga_mf <- CI_TCGA_mRNA %>% filter(Tumor_Sample_Barcode %in% both)

# Prepare gene expression matrix for PCA analysis
logfpkm_mf <- log2(mf_mrna_fpkm+1)
row.names(logfpkm_mf) <- row.names(mf_mrna_fpkm)

#Outliers 
gsg_mf <- goodSamplesGenes(t(logfpkm_mf))
summary(gsg_mf)
gsg_mf$allOK

table(gsg_mf$goodGenes)
table(gsg_mf$goodSamples)

# remove genes that are detectd as outliers
data_mf <- logfpkm_mf[gsg_mf$goodGenes == TRUE,]


# detect outlier samples - hierarchical clustering - method 1
htree_mf <- hclust(dist(t(data_mf)), method = "ward.D2")
plot(htree_mf)



#PCA
pca.res_mf <- pca(data_mf, metadata=ci_tcga_mf)

#Plot variance explained by each component
screeplot(pca.res_mf)

#Plot 2 selected components/eigenvectors
biplot(pca.res_mf,x="PC2", y="PC1")
biplot(pca.res_mf, hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
biplot(pca.res_mf, x="PC2", y="PC1",lab="", hline = 0, vline = 0,legendPosition = 'top') # Biplot without sample names
biplot(pca.res_mf, x="PC1", y="PC3",lab="",colby="TS", hline = 0, vline = 0,legendPosition = 'top') # Biplot with PC1 and PC3


#Plot several components
pairsplot(pca.res_mf)


# create dds
dds_mf <- DESeqDataSetFromMatrix(countData = round_data_mf,
                                 colData = col_metada_mf,
                                 design = ~ 1)

## remove all genes with low counts 
## suggested by WGCNA on RNAseq FAQ

dds75_mf <- dds_mf[rowSums(counts(dds_mf) >= 3) >= 24,]
nrow(dds75_mf)



# get normalized counts
norm.counts_mf <- assay(dds75_mf) %>% t()


# Network Construction 
# Choose a set of soft-thresholding powers
power_mf <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function

sft_mf <- pickSoftThreshold(norm.counts_mf,
                            powerVector = power_mf,
                            networkType = "signed",
                            verbose = 5)

sft.data_mf <- sft_mf$fitIndices

# visualization to pick power

a1_mf <- ggplot(sft.data_mf, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2_mf <- ggplot(sft.data_mf, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1_mf, a2_mf, nrow = 2)


# convert matrix to numeric
norm.counts_mf[] <- sapply(norm.counts_mf, as.numeric)

soft_power_mf <- 18
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet_mf <- blockwiseModules(norm.counts_mf,
                             maxBlockSize = 14000,
                             TOMType = "signed",
                             power = soft_power_mf,
                             mergeCutHeight = 0.25,
                             numericLabels = FALSE,
                             randomSeed = 1234,
                             verbose = 3)


cor <- temp_cor

# Module Eigengenes 
module_eigengenes_mf <- bwnet_mf$MEs


# Print out a preview
head(module_eigengenes_mf)


# get number of genes for each module
table(bwnet_mf$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_mf$dendrograms[[1]], cbind(bwnet_mf$unmergedColors, bwnet_mf$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


moduleColors_mf <- bwnet_mf$colors

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

# Relate modules to clinical traits
# module trait associations

# create traits file - binarize categorical variables
traits_mf <- col_metada_mf %>% 
  mutate(metastatisi = ifelse(grepl('YES', metastatic_disease_confirmed), 1, 0)) %>%
  mutate(gender = ifelse(grepl('FEMALE', gender), 1, 0)) %>%
  dplyr::select(1,3,6,23)

traits_mf <- traits_mf %>% mutate(status = ifelse(grepl('Alive', vital_status), 1, 0)) %>%
  mutate(otm = ifelse(grepl('Dead', history_other_malignancy), 1, 0))

traits_mf <- traits_mf %>% dplyr::select(1,4:6)


# binarize categorical variables


metastasis_site <- binarizeCategoricalColumns(col_metada_mf$metastasis_site,
                                              includePairwise = FALSE,
                                              includeLevelVsAll = TRUE,
                                              minCount = 1)

ethnicity <- binarizeCategoricalColumns(col_metada_mf$ethnicity,
                                   includePairwise = FALSE,
                                   includeLevelVsAll = TRUE,
                                   minCount = 1)

traits_mf <- cbind(traits_mf, metastasis_site,ethnicity)

#traits <- cbind(traits, severity.out)

# Define numbers of genes and samples
nSamples_mf <- nrow(norm.counts_mf)
nGenes_mf <- ncol(norm.counts_mf)


module.trait.corr_mf <- cor(module_eigengenes_mf, traits_mf, use = 'p')
module.trait.corr.pvals_mf <- corPvalueStudent(module.trait.corr_mf, nSamples_mf)



# visualize module-trait association as a heatmap

heatmap.data_mf <- merge(module_eigengenes_mf, traits_mf, by = 'row.names')

head(heatmap.data_mf)

heatmap.data_mf <- heatmap.data_mf[,-1]
colnames(heatmap.data_mf)[3:5] <- c("Gender","Metastasis","Status")
colnames(heatmap.data_mf)[7:12] <- c("Liver", "Lung", "Other", "Asian", "Black","White")

pdf("mf_wgcna.pdf")
CorLevelPlot(heatmap.data_mf,
             x = names(heatmap.data_mf)[c(3:5,7:12)],
             y = names(heatmap.data_mf)[1:2],
             col = c("blue1", "skyblue", "white", "pink", "red"),
             cexLabX = 0.6)

dev.off()

#See genes in module
module.gene.mapping <- as.data.frame(bwnet_mf$colors)
module.gene.mapping %>% 
  filter(`bwnet_mf$colors` == 'turquoise') %>% 
  rownames()

#Indentify hub genes
mf_s <- as.data.frame(traits_mf$status)
names(mf_s) <- "mf_s"

#names (colors) of the modules
modNames_mf <-  substring(names(module_eigengenes_mf), 3)

geneModuleMembership_mf <- as.data.frame(cor(norm.counts_mf, module_eigengenes_mf, use = "p"))
MMPvalue_mf <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_mf), nSamples_mf))

names(geneModuleMembership_mf) <- paste("MM", modNames_mf, sep="")
names(MMPvalue_mf) <- paste("p.MM", modNames_mf, sep="")

geneTraitSignificance_mf <- as.data.frame(cor(norm.counts_mf, mf_s, use = "p"))
GSPvalue_mf <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_mf), nSamples_mf))
GSPvalue_mf_adjust <- as.data.frame(p.adjust(as.matrix(GSPvalue_mf), method = "fdr"))

names(geneTraitSignificance_mf) <- paste("GS.", names(mf_s), sep="")
names(GSPvalue_mf) <- paste("p.GS.", names(mf_s), sep="")
names(GSPvalue_mf_adjust) <- paste("p.Adjust.GS.", names(mf_s), sep="")

module_mf <- "turquoise" #########################putting the color below the plot
column_mf <- match(module_mf, modNames_mf)
moduleGenes_mf <- moduleColors_mf==module_mf

verboseScatterplot(abs(geneModuleMembership_mf[moduleGenes_mf, column_mf]),
                   abs(geneTraitSignificance_mf[moduleGenes_mf, 1]),
                   xlab = paste("Module Membership in", module_mf, "module"),
                   ylab = "Gene significance for mf",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_mf)

geneInfo0_mf <- data.frame(EST = colnames(norm.counts_mf),
                          moduleColor = moduleColors_mf,
                          geneTraitSignificance_mf,
                          GSPvalue_mf, GSPvalue_mf_adjust)


modOrder_mf <- order(-abs(cor(module_eigengenes_mf, mf_s, use = "p")))
for (mod in 1:ncol(geneModuleMembership_mf))
{
  oldNames = names(geneInfo0_mf)
  geneInfo0_mf = data.frame(geneInfo0_mf, geneModuleMembership_mf[, modOrder_mf[mod]], 
                            MMPvalue_mf[, modOrder_mf[mod]]);
  names(geneInfo0_mf) = c(oldNames, paste("MM.", modNames_mf[modOrder_mf[mod]], sep=""),
                          paste("p.MM.", modNames_mf[modOrder_mf[mod]], sep=""))
}

geneOrder_mf <- order(geneInfo0_mf$moduleColor, -abs(geneInfo0_mf$GS.mf_s))
geneInfo_mf <- geneInfo0_mf[geneOrder_mf, ]

#Indentify Hub genes
turquoise_mf <- geneInfo_mf %>%
  filter(moduleColor=="turquoise") %>%
  filter(abs(GS.mf_s) >= 0.2) %>%
  filter(p.Adjust.GS.mf_s <= 0.05) %>%
  filter(abs(MM.turquoise) >= 0.8) %>%
  filter(p.MM.turquoise <= 0.05)


#####MPNST####
#Subbset for only ts
both <- intersect(rownames(CI_TCGA_mRNA),names(ts_geneexp[["MPNST"]]))

mpnst_mrna_fpkm <- mRNA_fpkm %>% 
  dplyr::select(all_of(both))

#Ci for ts
ci_tcga_mpnst <- CI_TCGA_mRNA %>% filter(Tumor_Sample_Barcode %in% both)

# Prepare gene expression matrix for PCA analysis: Step 1) get log TPMs (we add 1 TPM to each gene to avoid infinite values after log)
logfpkm_mpnst <- log2(mpnst_mrna_fpkm+1)
row.names(logfpkm_mpnst) <- row.names(mpnst_mrna_fpkm)

#Outliers 
gsg_mpnst <- goodSamplesGenes(t(logfpkm_mpnst))
summary(gsg_mpnst)
gsg_mpnst$allOK

table(gsg_mpnst$goodGenes)
table(gsg_mpnst$goodSamples)

# remove genes that are detectd as outliers
data_mpnst <- logfpkm_mpnst[gsg_mpnst$goodGenes == TRUE,]


# detect outlier samples - hierarchical clustering - method 1
htree_mpnst <- hclust(dist(t(data_mpnst)), method = "ward.D2")
plot(htree_mpnst)


#PCA
pca.res_mpnst <- pca(data_mpnst, metadata=ci_tcga_mpnst)

#Plot variance explained by each component
screeplot(pca.res_mpnst)

#Plot 2 selected components/eigenvectors
biplot(pca.res_mpnst)
biplot(pca.res_mpnst, hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
biplot(pca.res_mpnst, x="PC2", y="PC1",lab="", hline = 0, vline = 0,legendPosition = 'top') # Biplot without sample names
biplot(pca.res_mpnst, x="PC1", y="PC3",lab="",colby="TS", hline = 0, vline = 0,legendPosition = 'top') # Biplot with PC1 and PC3


#Plot several components
pairsplot(pca.res_mpnst)


#Prepare data for anlysis
round_data_mpnst <- round(data_mpnst)

col_metada_mpnst <- ci_tcga_mpnst[,c(7,8,10:18,32,34,46:51,67:69)]

# create dds
dds_mpnst <- DESeqDataSetFromMatrix(countData = round_data_mpnst,
                                 colData = col_metada_mpnst,
                                 design = ~ 1)

## remove all genes with counts low counts
## suggested by WGCNA on RNAseq FAQ

dds75_mpnst <- dds_mpnst[rowSums(counts(dds_mpnst) >= 2) >= 6,]
nrow(dds75_mpnst)


# perform variance stabilization
dds_norm_mpnst <- varianceStabilizingTransformation(dds75_mpnst)

# get normalized counts
norm.counts_mpnst <- assay(dds_norm_mpnst) %>% t()


# Network Construction 
# Choose a set of soft-thresholding powers
power_mpnst <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function

sft_mpnst <- pickSoftThreshold(norm.counts_mpnst,
                            powerVector = power_mpnst,
                            networkType = "signed",
                            verbose = 5)

sft.data_mpnst <- sft_mpnst$fitIndices

# visualization to pick power

a1_mpnst <- ggplot(sft.data_mpnst, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2_mpnst <- ggplot(sft.data_mpnst, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1_mpnst, a2_mpnst, nrow = 2)


# convert matrix to numeric
norm.counts_mpnst[] <- sapply(norm.counts_mpnst, as.numeric)

soft_power_mpnst <- 22
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet_mpnst <- blockwiseModules(norm.counts_mpnst,
                             maxBlockSize = 14000,
                             TOMType = "signed",
                             power = soft_power_mpnst,
                             mergeCutHeight = 0.25,
                             numericLabels = FALSE,
                             randomSeed = 1234,
                             verbose = 3)


cor <- temp_cor

# Module Eigengenes 
module_eigengenes_mpnst <- bwnet_mpnst$MEs


# Print out a preview
head(module_eigengenes_mpnst)


# get number of genes for each module
table(bwnet_mpnst$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_mpnst$dendrograms[[1]], cbind(bwnet_mpnst$unmergedColors, bwnet_mpnst$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


moduleColors_mpnst <- bwnet_mpnst$colors

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

#Relate modules to clinical traits
# module trait associations



# create traits file - binarize categorical variables
traits_mpnst <- col_metada_mpnst %>% 
  mutate(metastatisi = ifelse(grepl('YES', metastatic_disease_confirmed), 1, 0)) %>%
  mutate(gender = ifelse(grepl('FEMALE', gender), 1, 0)) %>%
  dplyr::select(1,3,6,23)

traits_mpnst <- traits_mpnst %>% mutate(status = ifelse(grepl('Alive', vital_status), 1, 0)) %>%
  mutate(otm = ifelse(grepl('Dead', history_other_malignancy), 1, 0))

traits_mpnst <- traits_mpnst %>% dplyr::select(1,4:6)


# binarize categorical variables


metastasis_site <- binarizeCategoricalColumns(col_metada_mpnst$metastasis_site,
                                              includePairwise = FALSE,
                                              includeLevelVsAll = TRUE,
                                              minCount = 1)

ethnicity <- binarizeCategoricalColumns(col_metada_mpnst$ethnicity,
                                   includePairwise = FALSE,
                                   includeLevelVsAll = TRUE,
                                   minCount = 1)

traits_mpnst <- cbind(traits_mpnst, metastasis_site, ethnicity)


# Define numbers of genes and samples
nSamples_mpnst <- nrow(norm.counts_mpnst)
nGenes_mpnst <- ncol(norm.counts_mpnst)


module.trait.corr_mpnst <- cor(module_eigengenes_mpnst, traits_mpnst, use = 'p')
module.trait.corr.pvals_mpnst <- corPvalueStudent(module.trait.corr_mpnst, nSamples_mpnst)



# visualize module-trait association as a heatmap

heatmap.data_mpnst <- merge(module_eigengenes_mpnst, traits_mpnst, by = 'row.names')

head(heatmap.data_mpnst)

heatmap.data_mpnst <- heatmap.data_mpnst[,-1]
colnames(heatmap.data_mpnst)[90:92] <- c("Gender","Metastasis","Status")
colnames(heatmap.data_mpnst)[95] <- c("White")

pdf("mpnst_wgcna.pdf")
CorLevelPlot(heatmap.data_mpnst,
             x = names(heatmap.data_mpnst)[c(90:92,95)],
             y = names(heatmap.data_mpnst)[1:89], cexLabX = 0.8,
             col = c("blue1", "skyblue", "white", "pink", "red"),
             cexLabY = 0.3, cexCorval = 0.4
             )

dev.off()

#Visualize the genes in a module
module.gene.mapping <- as.data.frame(bwnet_mpnst$colors)
module.gene.mapping %>% 
  filter(`bwnet_mpnst$colors` == 'salmon4') %>% 
  rownames()

##Indetify Hub Genes
#Metastatis as clinical trait
mpnst_m <- as.data.frame(traits_mpnst$Metastasis)
names(mpnst_m) <- "mpnst_m"

#names (colors) of the modules
modNames_mpnst <-  substring(names(module_eigengenes_mpnst), 3)

geneModuleMembership_mpnst <- as.data.frame(cor(norm.counts_mpnst, module_eigengenes_mpnst, use = "p"))
MMPvalue_mpnst <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_mpnst), nSamples_mpnst))


names(geneModuleMembership_mpnst) <- paste("MM", modNames_mpnst, sep="")
names(MMPvalue_mpnst) <- paste("p.MM", modNames_mpnst, sep="")

geneTraitSignificance_mpnst <- as.data.frame(cor(norm.counts_mpnst, mpnst_m, use = "p"))
GSPvalue_mpnst <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_mpnst), nSamples_mpnst))
GSPvalue_mpnst_adjust <- as.data.frame(p.adjust(as.matrix(GSPvalue_mpnst), method="fdr"))

names(geneTraitSignificance_mpnst) <- paste("GS.", names(mpnst_m), sep="")
names(GSPvalue_mpnst) <- paste("p.GS.", names(mpnst_m), sep="")
names(GSPvalue_mpnst_adjust) <- paste("p.Adjust.GS.", names(mpnst_m), sep="")

module_mpnst <- "salmon4" #########################putting the color below the plot
column_mpnst <- match(module_mpnst, modNames_mpnst)
moduleGenes_mpnst <- moduleColors_mpnst==module_mpnst

verboseScatterplot(abs(geneModuleMembership_mpnst[moduleGenes_mpnst, column_mpnst]),
                   abs(geneTraitSignificance_mpnst[moduleGenes_mpnst, 1]),
                   xlab = paste("Module Membership in", module_mpnst, "module"),
                   ylab = "Gene significance for mpnst",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

geneInfo0_mpnst <- data.frame(EST = colnames(norm.counts_mpnst),
                             moduleColor = moduleColors_mpnst,
                             geneTraitSignificance_mpnst,
                             GSPvalue_mpnst,
                             GSPvalue_mpnst_adjust)


modOrder_mpnst <- order(-abs(cor(module_eigengenes_mpnst, mpnst_m, use = "p")))
for (mod in 1:ncol(geneModuleMembership_mpnst))
{
  oldNames = names(geneInfo0_mpnst)
  geneInfo0_mpnst = data.frame(geneInfo0_mpnst, geneModuleMembership_mpnst[, modOrder_mpnst[mod]], 
                               MMPvalue_mpnst[, modOrder_mpnst[mod]]);
  names(geneInfo0_mpnst) = c(oldNames, paste("MM.", modNames_mpnst[modOrder_mpnst[mod]], sep=""),
                             paste("p.MM.", modNames_mpnst[modOrder_mpnst[mod]], sep=""))
}
geneOrder_mpnst <- order(geneInfo0_mpnst$moduleColor, -abs(geneInfo0_mpnst$GS.mpnst_m))
geneInfo_mpnst <- geneInfo0_mpnst[geneOrder_mpnst, ]

#Indentify Hub genes
salmon4_mpnst <- geneInfo_mpnst %>%
  filter(moduleColor=="salmon4") %>%
  filter(abs(GS.mpnst_m) >= 0.2) %>%
  filter(p.Adjust.GS.mpnst_m <= 0.05) %>%
  filter(abs(MM.salmon4) >= 0.8) %>%
  filter(p.MM.salmon4 <= 0.05)

#Status
mpnst_s <- as.data.frame(traits_mpnst$Status)
names(mpnst_s) <- "mpnst_s"

#names (colors) of the modules
modNames_mpnst_s <-  substring(names(module_eigengenes_mpnst), 3)

geneModuleMembership_mpnst_s <- as.data.frame(cor(norm.counts_mpnst, module_eigengenes_mpnst, use = "p"))
MMPvalue_mpnst_s <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_mpnst_s), nSamples_mpnst))

names(geneModuleMembership_mpnst_s) <- paste("MM", modNames_mpnst_s, sep="")
names(MMPvalue_mpnst_s) <- paste("p.MM", modNames_mpnst_s, sep="")

geneTraitSignificance_mpnst_s <- as.data.frame(cor(norm.counts_mpnst, mpnst_s, use = "p"))
GSPvalue_mpnst_s <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_mpnst_s), nSamples_mpnst))
GSPvalue_mpnst_adjust_s <- as.data.frame(p.adjust(as.matrix(GSPvalue_mpnst_s), method="fdr"))

names(geneTraitSignificance_mpnst_s) <- paste("GS.", names(mpnst_s), sep="")
names(GSPvalue_mpnst_adjust_s) <- paste("p.GS.", names(mpnst_s), sep="")
names(GSPvalue_mpnst_adjust_s) <- paste("p.Adjust.GS.", names(mpnst_s), sep="")

module_mpnst <- "black" #########################putting the color below the plot
column_mpnst <- match(module_mpnst, modNames_mpnst)
moduleGenes_mpnst <- moduleColors_mpnst==module_mpnst

#Plot Module Membership Vs Gene significance
verboseScatterplot(abs(geneModuleMembership_mpnst_s[moduleGenes_mpnst, column_mpnst]),
                   abs(geneTraitSignificance_mpnst_s[moduleGenes_mpnst, 1]),
                   xlab = paste("Module Membership in", module_mpnst, "module"),
                   ylab = "Gene significance for mpnst",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_mpnst)

#Create dataframe with all the info
geneInfo0_mpnst_s <- data.frame(EST = colnames(norm.counts_mpnst),
                               moduleColor = moduleColors_mpnst,
                               geneTraitSignificance_mpnst_s,
                               GSPvalue_mpnst_s,
                               GSPvalue_mpnst_adjust_s)


modOrder_mpnst_s <- order(-abs(cor(module_eigengenes_mpnst, mpnst_s, use = "p")))
for (mod in 1:ncol(geneModuleMembership_mpnst_s)){
  oldNames = names(geneInfo0_mpnst_s)
  geneInfo0_mpnst_s = data.frame(geneInfo0_mpnst_s, geneModuleMembership_mpnst_s[, modOrder_mpnst_s[mod]], 
                                 MMPvalue_mpnst_s[, modOrder_mpnst_s[mod]]);
  names(geneInfo0_mpnst_s) = c(oldNames, paste("MM.", modNames_mpnst_s[modOrder_mpnst_s[mod]], sep=""),
                               paste("p.MM.", modNames_mpnst_s[modOrder_mpnst_s[mod]], sep=""))
}
geneOrder_mpnst_s <- order(geneInfo0_mpnst_s$moduleColor, -abs(geneInfo0_mpnst_s$GS.mpnst_s))
geneInfo_mpnst_s <- geneInfo0_mpnst_s[geneOrder_mpnst_s, ]

#Indentify Hub genes
black_mpnst <- geneInfo_mpnst_s %>%
  filter(moduleColor=="black") %>%
  filter(abs(GS.mpnst_s) >= 0.2) %>%
  filter(p.Adjust.GS.mpnst_s <= 0.05) %>%
  filter(abs(MM.black) >= 0.8) %>%
  filter(p.MM.black <= 0.05)

#Gender
mpnst_g <- as.data.frame(traits_mpnst$Gender)
names(mpnst_g) <- "mpnst_g"

#names (colors) of the modules
modNames_mpnst_g <-  substring(names(module_eigengenes_mpnst), 3)

geneModuleMembership_mpnst_g <- as.data.frame(cor(norm.counts_mpnst, module_eigengenes_mpnst, use = "p"))
MMPvalue_mpnst_g <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_mpnst_g), nSamples_mpnst))

names(geneModuleMembership_mpnst_g) <- paste("MM", modNames_mpnst_g, sep="")
names(MMPvalue_mpnst_g) <- paste("p.MM", modNames_mpnst_g, sep="")

geneTraitSignificance_mpnst_g <- as.data.frame(cor(norm.counts_mpnst, mpnst_g, use = "p"))
GSPvalue_mpnst_g <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_mpnst_g), nSamples_mpnst))
GSPvalue_mpnst_adjust_g <- as.data.frame(p.adjust(as.matrix(GSPvalue_mpnst_g), method="fdr"))

names(geneTraitSignificance_mpnst_g) <- paste("GS.", names(mpnst_g), sep="")
names(GSPvalue_mpnst_adjust_g) <- paste("p.GS.", names(mpnst_g), sep="")
names(GSPvalue_mpnst_adjust_g) <- paste("p.Adjust.GS.", names(mpnst_g), sep="")

module_mpnst <- "midnightblue" #########################putting the color below the plot
column_mpnst <- match(module_mpnst, modNames_mpnst)
moduleGenes_mpnst <- moduleColors_mpnst==module_mpnst

#Plot Module Membership Vs Gene significance
verboseScatterplot(abs(geneModuleMembership_mpnst_g[moduleGenes_mpnst, column_mpnst]),
                   abs(geneTraitSignificance_mpnst_g[moduleGenes_mpnst, 1]),
                   xlab = paste("Module Membership in", module_mpnst, "module"),
                   ylab = "Gene significance for mpnst",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_mpnst)

#Create dataframe with all the info
geneInfo0_mpnst_g <- data.frame(EST = colnames(norm.counts_mpnst),
                                moduleColor = moduleColors_mpnst,
                                geneTraitSignificance_mpnst_g,
                                GSPvalue_mpnst_g,
                                GSPvalue_mpnst_adjust_s)


modOrder_mpnst_g <- order(-abs(cor(module_eigengenes_mpnst, mpnst_g, use = "p")))
for (mod in 1:ncol(geneModuleMembership_mpnst_g)){
  oldNames = names(geneInfo0_mpnst_g)
  geneInfo0_mpnst_g = data.frame(geneInfo0_mpnst_g, geneModuleMembership_mpnst_g[, modOrder_mpnst_g[mod]], 
                                 MMPvalue_mpnst_g[, modOrder_mpnst_g[mod]]);
  names(geneInfo0_mpnst_g) = c(oldNames, paste("MM.", modNames_mpnst_g[modOrder_mpnst_g[mod]], sep=""),
                               paste("p.MM.", modNames_mpnst_g[modOrder_mpnst_g[mod]], sep=""))
}
geneOrder_mpnst_g <- order(geneInfo0_mpnst_g$moduleColor, -abs(geneInfo0_mpnst_g$GS.mpnst_g))
geneInfo_mpnst_g <- geneInfo0_mpnst_g[geneOrder_mpnst_g, ]

#Indentify Hub genes
mb_mpnst <- geneInfo_mpnst_g %>%
  filter(moduleColor=="midnightblue") %>%
  filter(abs(GS.mpnst_g) >= 0.2) %>%
  filter(p.Adjust.GS.mpnst_g <= 0.05) %>%
  filter(abs(MM.midnightblue) >= 0.8) %>%
  filter(p.MM.midnightblue <= 0.05)

#Ethinicty
mpnst_e <- as.data.frame(traits_mpnst$White)
names(mpnst_e) <- "mpnst_e"

#names (colors) of the modules
modNames_mpnst_e <-  substring(names(module_eigengenes_mpnst), 3)

geneModuleMembership_mpnst_e <- as.data.frame(cor(norm.counts_mpnst, module_eigengenes_mpnst, use = "p"))
MMPvalue_mpnst_e <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_mpnst_e), nSamples_mpnst))

names(geneModuleMembership_mpnst_e) <- paste("MM", modNames_mpnst_e, sep="")
names(MMPvalue_mpnst_e) <- paste("p.MM", modNames_mpnst_e, sep="")

geneTraitSignificance_mpnst_e <- as.data.frame(cor(norm.counts_mpnst, mpnst_e, use = "p"))
GSPvalue_mpnst_e <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_mpnst_e), nSamples_mpnst))
GSPvalue_mpnst_adjust_e <- as.data.frame(p.adjust(as.matrix(GSPvalue_mpnst_e), method="fdr"))

names(geneTraitSignificance_mpnst_e) <- paste("GS.", names(mpnst_e), sep="")
names(GSPvalue_mpnst_adjust_e) <- paste("p.GS.", names(mpnst_e), sep="")
names(GSPvalue_mpnst_adjust_e) <- paste("p.Adjust.GS.", names(mpnst_e), sep="")

module_mpnst <- "darkgreen" #########################putting the color below the plot
column_mpnst <- match(module_mpnst, modNames_mpnst)
moduleGenes_mpnst <- moduleColors_mpnst==module_mpnst

#Plot Module Membership Vs Gene significance
verboseScatterplot(abs(geneModuleMembership_mpnst_e[moduleGenes_mpnst, column_mpnst]),
                   abs(geneTraitSignificance_mpnst_e[moduleGenes_mpnst, 1]),
                   xlab = paste("Module Membership in", module_mpnst, "module"),
                   ylab = "Gene significance for mpnst",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_mpnst)

#Create dataframe with all the info
geneInfo0_mpnst_e <- data.frame(EST = colnames(norm.counts_mpnst),
                                moduleColor = moduleColors_mpnst,
                                geneTraitSignificance_mpnst_e,
                                GSPvalue_mpnst_e,
                                GSPvalue_mpnst_adjust_s)


modOrder_mpnst_e <- order(-abs(cor(module_eigengenes_mpnst, mpnst_e, use = "p")))
for (mod in 1:ncol(geneModuleMembership_mpnst_e)){
  oldNames = names(geneInfo0_mpnst_e)
  geneInfo0_mpnst_e = data.frame(geneInfo0_mpnst_e, geneModuleMembership_mpnst_e[, modOrder_mpnst_e[mod]], 
                                 MMPvalue_mpnst_e[, modOrder_mpnst_e[mod]]);
  names(geneInfo0_mpnst_e) = c(oldNames, paste("MM.", modNames_mpnst_e[modOrder_mpnst_e[mod]], sep=""),
                               paste("p.MM.", modNames_mpnst_e[modOrder_mpnst_e[mod]], sep=""))
}
geneOrder_mpnst_e <- order(geneInfo0_mpnst_e$moduleColor, -abs(geneInfo0_mpnst_e$GS.mpnst_e))
geneInfo_mpnst_e <- geneInfo0_mpnst_e[geneOrder_mpnst_e, ]

#Indentify Hub genes
dg_mpnst <- geneInfo_mpnst_e %>%
  filter(moduleColor=="darkgreen") %>%
  filter(abs(GS.mpnst_e) >= 0.2) %>%
  filter(p.Adjust.GS.mpnst_e <= 0.05) %>%
  filter(abs(MM.darkgreen) >= 0.8) %>%
  filter(p.MM.darkgreen <= 0.05)

maroon_mpnst <- geneInfo_mpnst_e %>%
  filter(moduleColor=="maroon") %>%
  filter(abs(GS.mpnst_e) >= 0.2) %>%
  filter(p.Adjust.GS.mpnst_e <= 0.05) %>%
  filter(abs(MM.maroon) >= 0.8) %>%
  filter(p.MM.maroon <= 0.05)

#SAve hub genes data
write.xlsx(maroon_mpnst, "mpnst_wgna_m.xlsx")
write.xlsx(dg_mpnst, "mpnst_wgna_dg.xlsx")
write.xlsx(mb_mpnst, "mpnst_wgna_mb.xlsx")
write.xlsx(black_mpnst, "mpnst_wgna_b.xlsx")
write.xlsx(salmon4_mpnst, "mpnst_wgna_salmon.xlsx")


#####UPS####
#Subbset for only ts
both <- intersect(rownames(CI_TCGA_mRNA),names(ts_geneexp[["UPS"]]))

ups_mrna_fpkm <- mRNA_fpkm %>% 
  dplyr::select(all_of(both))

#Ci for ts
ci_tcga_ups <- CI_TCGA_mRNA %>% filter(Tumor_Sample_Barcode %in% both)

# Prepare gene expression matrix for PCA analysis
logfpkm_ups <- log2(ups_mrna_fpkm+1)
row.names(logfpkm_ups) <- row.names(ups_mrna_fpkm)

#Outliers 
gsg_ups <- goodSamplesGenes(t(logfpkm_ups))
summary(gsg_ups)
gsg_ups$allOK

table(gsg_ups$goodGenes)
table(gsg_ups$goodSamples)

# remove genes that are detectd as outliers
data_ups <- logfpkm_ups[gsg_ups$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree_ups <- hclust(dist(t(data_ups)), method = "ward.D2")
plot(htree_ups)



#PCA
pca.res_ups <- pca(data_ups, metadata=ci_tcga_ups)

#Plot variance explained by each component
screeplot(pca.res_ups)

#Plot 2 selected components/eigenvectors
biplot(pca.res_ups)
biplot(pca.res_ups, hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
biplot(pca.res_ups, x="PC2", y="PC1",lab="", hline = 0, vline = 0,legendPosition = 'top') # Biplot without sample names
biplot(pca.res_ups, x="PC1", y="PC3",lab="",colby="TS", hline = 0, vline = 0,legendPosition = 'top') # Biplot with PC1 and PC3


#Plot several components
pairsplot(pca.res_ups)

#Prepare Data for analysis
round_data_ups <- round(data_ups)

col_metada_ups <- ci_tcga_ups[,c(7,8,10:18,32,34,46:51,67:69)]

# create dds
dds_ups <- DESeqDataSetFromMatrix(countData = round_data_ups,
                                    colData = col_metada_ups,
                                    design = ~ 1)

## remove all genes with low counts 
## suggested by WGCNA on RNAseq FAQ

dds75_ups <- dds_ups[rowSums(counts(dds_ups) >= 3) >= 16,]
nrow(dds75_ups)


# perform variance stabilization
dds_norm_ups <- varianceStabilizingTransformation(dds75_ups)

# get normalized counts
norm.counts_ups <- assay(dds_norm_ups) %>% t()


#Network Construction 
# Choose a set of soft-thresholding powers
power_ups <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function

sft_ups <- pickSoftThreshold(norm.counts_ups,
                               powerVector = power_ups,
                               networkType = "signed",
                               verbose = 5)

sft.data_ups <- sft_ups$fitIndices

# visualization to pick power

a1_ups <- ggplot(sft.data_ups, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2_ups <- ggplot(sft.data_ups, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1_ups, a2_ups, nrow = 2)


# convert matrix to numeric
norm.counts_ups[] <- sapply(norm.counts_ups, as.numeric)

soft_power_ups <- 12
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet_ups <- blockwiseModules(norm.counts_ups,
                                maxBlockSize = 14000,
                                TOMType = "signed",
                                power = soft_power_ups,
                                mergeCutHeight = 0.25,
                                numericLabels = FALSE,
                                randomSeed = 1234,
                                verbose = 3)


cor <- temp_cor

#Module Eigengenes 
module_eigengenes_ups <- bwnet_ups$MEs


# Print out a preview
head(module_eigengenes_ups)


# get number of genes for each module
table(bwnet_ups$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_ups$dendrograms[[1]], cbind(bwnet_ups$unmergedColors, bwnet_ups$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


moduleColors_ups <- bwnet_ups$colors

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

# Relate modules to traits
# module trait associations



# create traits file - binarize categorical variables
traits_ups <- col_metada_ups %>% 
  mutate(metastatisi = ifelse(grepl('YES', metastatic_disease_confirmed), 1, 0)) %>%
  mutate(gender = ifelse(grepl('FEMALE', gender), 1, 0)) %>%
  dplyr::select(1,3,6,23)

traits_ups <- traits_ups %>% mutate(status = ifelse(grepl('Alive', vital_status), 1, 0)) %>%
  mutate(otm = ifelse(grepl('Dead', history_other_malignancy), 1, 0))

traits_ups <- traits_ups %>% dplyr::select(1,4:6)


# binarize categorical variables


metastasis_site <- binarizeCategoricalColumns(col_metada_ups$metastasis_site,
                                              includePairwise = FALSE,
                                              includeLevelVsAll = TRUE,
                                              minCount = 1)

ethnicity<- binarizeCategoricalColumns(col_metada_ups$ethnicity,
                                   includePairwise = FALSE,
                                   includeLevelVsAll = TRUE,
                                   minCount = 1)

traits_ups <- cbind(traits_ups, metastasis_site, ethnicity)


# Define numbers of genes and samples
nSamples_ups <- nrow(norm.counts_ups)
nGenes_ups <- ncol(norm.counts_ups)


module.trait.corr_ups <- cor(module_eigengenes_ups, traits_ups, use = 'p')
module.trait.corr.pvals_ups <- corPvalueStudent(module.trait.corr_ups, nSamples_ups)



# visualize module-trait association as a heatmap

heatmap.data_ups <- merge(module_eigengenes_ups, traits_ups, by = 'row.names')

head(heatmap.data_ups)

heatmap.data_ups <- heatmap.data_ups[,-1]
colnames(heatmap.data_ups)[9:11] <- c("Gender","Metastasis","Status")
colnames(heatmap.data_ups)[13:16] <- c("Lung","Other","Black","White") 

pdf("ups_wgcna.pdf")
CorLevelPlot(heatmap.data_ups,
             x = names(heatmap.data_ups)[c(9:11,13:16)],
             y = names(heatmap.data_ups)[1:8], cexLabX = 0.8,
             col = c("blue1", "skyblue", "white", "pink", "red")
             )

dev.off()

#See genes in the module
module.gene.mapping <- as.data.frame(bwnet_ups$colors)
module.gene.mapping %>% 
  filter(`bwnet_ups$colors` == 'turquoise') %>% 
  rownames()

#Indentify Hub genes
#Metastasis
ups_m <- as.data.frame(traits_ups$Metastasis)
names(ups_m) <- "ups_m"

#names (colors) of the modules
modNames_ups <-  substring(names(module_eigengenes_ups), 3)

geneModuleMembership_ups <- as.data.frame(cor(norm.counts_ups, module_eigengenes_ups, use = "p"))
MMPvalue_ups <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_ups), nSamples_ups))

names(geneModuleMembership_ups) <- paste("MM", modNames_ups, sep="")
names(MMPvalue_ups) <- paste("p.MM", modNames_ups, sep="")

geneTraitSignificance_ups <- as.data.frame(cor(norm.counts_ups, ups_m, use = "p"))
GSPvalue_ups <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_ups), nSamples_ups))
GSPvalue_ups_adjust <- as.data.frame(p.adjust(as.matrix(GSPvalue_ups), method = "fdr"))

names(geneTraitSignificance_ups) <- paste("GS.", names(ups_m), sep="")
names(GSPvalue_ups) <- paste("p.GS.", names(ups_m), sep="")
names(GSPvalue_ups_adjust) <- paste("p.Adjust.GS.", names(ups_m), sep="")

module_ups <- "yellow" #########################putting the color below the plot
column_ups <- match(module_ups, modNames_ups)
moduleGenes_ups <- moduleColors_ups==module_ups

verboseScatterplot(abs(geneModuleMembership_ups[moduleGenes_ups, column_ups]),
                   abs(geneTraitSignificance_ups[moduleGenes_ups, 1]),
                   xlab = paste("Module Membership in", module_ups, "module"),
                   ylab = "Gene significance for ups",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_ups)

geneInfo0_ups <- data.frame(EST = colnames(norm.counts_ups),
                           moduleColor = moduleColors_ups,
                           geneTraitSignificance_ups,
                           GSPvalue_ups, 
                           GSPvalue_ups_adjust)


modOrder_ups <- order(-abs(cor(module_eigengenes_ups, ups_m, use = "p")))
for (mod in 1:ncol(geneModuleMembership_ups))
{
  oldNames = names(geneInfo0_ups)
  geneInfo0_ups = data.frame(geneInfo0_ups, geneModuleMembership_ups[, modOrder_ups[mod]], 
                             MMPvalue_ups[, modOrder_ups[mod]]);
  names(geneInfo0_ups) = c(oldNames, paste("MM.", modNames_ups[modOrder_ups[mod]], sep=""),
                           paste("p.MM.", modNames_ups[modOrder_ups[mod]], sep=""))
}
geneOrder_ups <- order(geneInfo0_ups$moduleColor, -abs(geneInfo0_ups$GS.ups_m))
geneInfo_ups <- geneInfo0_ups[geneOrder_ups, ]

yellow_ups <- geneInfo_ups %>%
  filter(moduleColor=="yellow") %>%
  filter(abs(GS.ups_m) >= 0.2) %>%
  filter(p.Adjust.GS.ups_m <= 0.05) %>%
  filter(abs(MM.yellow) >= 0.8) %>%
  filter(p.MM.yellow <= 0.05)

green_ups <- geneInfo_ups %>%
  filter(moduleColor=="green") %>%
  filter(abs(GS.ups_m) >= 0.2) %>%
  filter(p.Adjust.GS.ups_m <= 0.05) %>%
  filter(abs(MM.green) >= 0.8) %>%
  filter(p.MM.green <= 0.05)

#Status
ups_s <- as.data.frame(traits_ups$Status)
names(ups_s) <- "ups_s"

#names (colors) of the modules
modNames_ups <-  substring(names(module_eigengenes_ups), 3)

geneModuleMembership_ups <- as.data.frame(cor(norm.counts_ups, module_eigengenes_ups, use = "p"))
MMPvalue_ups <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_ups), nSamples_ups))
names(geneModuleMembership_ups) <- paste("MM", modNames_ups, sep="")
names(MMPvalue_ups) <- paste("p.MM", modNames_ups, sep="")

geneTraitSignificance_ups_s <- as.data.frame(cor(norm.counts_ups, ups_s, use = "p"))
GSPvalue_ups_s <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_ups_s), nSamples_ups))
GSPvalue_ups_adjust_s <- as.data.frame(p.adjust(as.matrix(GSPvalue_ups_s), method = "fdr"))

names(geneTraitSignificance_ups_s) <- paste("GS.", names(ups_s), sep="")
names(GSPvalue_ups_s) <- paste("p.GS.", names(ups_s), sep="")
names(GSPvalue_ups_adjust_s) <- paste("p.Adjust.GS.", names(ups_s), sep="")

module_ups <- "turquoise" #########################putting the color below the plot
column_ups <- match(module_ups, modNames_ups)
moduleGenes_ups <- moduleColors_ups==module_ups

verboseScatterplot(abs(geneModuleMembership_ups[moduleGenes_ups, column_ups]),
                   abs(geneTraitSignificance_ups[moduleGenes_ups, 1]),
                   xlab = paste("Module Membership in", module_ups, "module"),
                   ylab = "Gene significance for ups",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_ups)

geneInfo0_ups <- data.frame(EST = colnames(norm.counts_ups),
                           moduleColor = moduleColors_ups,
                           geneTraitSignificance_ups,
                           GSPvalue_ups, 
                           GSPvalue_ups_adjust)


modOrder_ups <- order(-abs(cor(module_eigengenes_ups, ups_s, use = "p")))
for (mod in 1:ncol(geneModuleMembership_ups))
{
  oldNames = names(geneInfo0_ups)
  geneInfo0_ups = data.frame(geneInfo0_ups, geneModuleMembership_ups[, modOrder_ups[mod]], 
                             MMPvalue_ups[, modOrder_ups[mod]]);
  names(geneInfo0_ups) = c(oldNames, paste("MM.", modNames_ups[modOrder_ups[mod]], sep=""),
                           paste("p.MM.", modNames_ups[modOrder_ups[mod]], sep=""))
}
geneOrder_ups <- order(geneInfo0_ups$moduleColor, -abs(geneInfo0_ups$GS.ups_s))
geneInfo_ups <- geneInfo0_ups[geneOrder_ups, ]

turquoise_ups <- geneInfo_ups %>%
  filter(moduleColor=="turquoise") %>%
  filter(abs(GS.ups_s) >= 0.2) %>%
  filter(p.Adjust.GS.ups_s <= 0.05) %>%
  filter(abs(MM.turquoise) >= 0.8) %>%
  filter(p.MM.turquoise <= 0.05)

black_ups <- geneInfo_ups %>%
  filter(moduleColor=="black") %>%
  filter(abs(GS.ups_s) >= 0.2) %>%
  filter(p.Adjust.GS.ups_s <= 0.05) %>%
  filter(abs(MM.black) >= 0.8) %>%
  filter(p.MM.black <= 0.05)

#####Synovial Sarcoma####
#Subbset for only ts
both <- intersect(rownames(CI_TCGA_mRNA),names(ts_geneexp[["SS"]]))

ss_mrna_fpkm <- mRNA_fpkm %>% 
  dplyr::select(all_of(both))

#Ci for ts
ci_tcga_ss <- CI_TCGA_mRNA %>% filter(Tumor_Sample_Barcode %in% both)

# Prepare gene expression matrix for PCA analysis: Step 1) get log TPMs (we add 1 TPM to each gene to avoid infinite values after log)
logfpkm_ss <- log2(ss_mrna_fpkm+1)
row.names(logfpkm_ss) <- row.names(ss_mrna_fpkm)

#Outliers 
gsg_ss <- goodSamplesGenes(t(logfpkm_ss))
summary(gsg_ss)
gsg_ss$allOK

table(gsg_ss$goodGenes)
table(gsg_ss$goodSamples)

# remove genes that are detectd as outliers
data_ss <- logfpkm_ss[gsg_ss$goodGenes == TRUE,]
#keep <- filterByExpr(data_mRNA_logfpkm, group=condition_ts_mRNA)
#data_mRNA_logfpkm <- data_mRNA_logfpkm[keep,]

# detect outlier samples - hierarchical clustering - method 1
htree_ss <- hclust(dist(t(data_ss)), method = "ward.D2")
plot(htree_ss)


data_ss <- data_ss %>% dplyr::select(!(c("TCGA-DX-A7EQ", "TCGA-Z4-AAPF")))
ci_tcga_ss <- ci_tcga_ss %>% filter(!(Tumor_Sample_Barcode %in% c("TCGA-DX-A7EQ", "TCGA-Z4-AAPF")))
#PCA
pca.res_ss <- pca(data_ss, metadata=ci_tcga_ss)

#Plot variance explained by each component
screeplot(pca.res_ss)

#Plot 2 selected components/eigenvectors
biplot(pca.res_ss)
biplot(pca.res_ss, hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
biplot(pca.res_ss, x="PC2", y="PC1",lab="", hline = 0, vline = 0,legendPosition = 'top') # Biplot without sample names
biplot(pca.res_ss, x="PC1", y="PC3",lab="",colby="TS", hline = 0, vline = 0,legendPosition = 'top') # Biplot with PC1 and PC3


#Plot several components
pairsplot(pca.res_ss)


PC_genes_ss <- pca.res_ss$loadings
PC1_genes_ss <- PC_genes[order(PC_genes_ss$PC1, decreasing=T),]
head(PC1_genes_ss)
tail(PC1_genes_ss)
plotloadings(pca.res_ss, components = c("PC1", "PC2", "PC3"),rangeRetain =0.1) #retaining 1% of the loadings per PC

round_data_ss <- round(data_ss)

col_metada_ss <- ci_tcga_ss[,c(7,8,10:18,32,34,46:51,67:69)]

# create dds
dds_ss <- DESeqDataSetFromMatrix(countData = round_data_ss,
                                    colData = col_metada_ss,
                                    design = ~ 1)

## remove all genes with low counts

dds75_ss <- dds_ss[rowSums(counts(dds_ss) >= 2) >= 6,]
nrow(dds75_ss)


# perform variance stabilization
dds_norm_ss <- varianceStabilizingTransformation(dds75_ss)

# get normalized counts
norm.counts_ss <- assay(dds_norm_ss) %>% t()


# Network Construction 
# Choose a set of soft-thresholding powers
power_ss <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function

sft_ss <- pickSoftThreshold(norm.counts_ss,
                               powerVector = power_ss,
                               networkType = "signed",
                               verbose = 5)

sft.data_ss <- sft_ss$fitIndices

# visualization to pick power

a1_ss <- ggplot(sft.data_ss, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2_ss <- ggplot(sft.data_ss, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1_ss, a2_ss, nrow = 2)


# convert matrix to numeric
norm.counts_ss[] <- sapply(norm.counts_ss, as.numeric)

soft_power_ss <- 12
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet_ss <- blockwiseModules(norm.counts_ss,
                                maxBlockSize = 14000,
                                TOMType = "signed",
                                power = soft_power_ss,
                                mergeCutHeight = 0.25,
                                numericLabels = FALSE,
                                randomSeed = 1234,
                                verbose = 3)


cor <- temp_cor

#Module Eigengenes 
module_eigengenes_ss <- bwnet_ss$MEs


# Print out a preview
head(module_eigengenes_ss)


# get number of genes for each module
table(bwnet_ss$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_ss$dendrograms[[1]], cbind(bwnet_ss$unmergedColors, bwnet_ss$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


moduleColors_ss <- bwnet_ss$colors

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

# Relate modules to clinical traits
# module trait associations



# create traits file - binarize categorical variables
traits_ss <- col_metada_ss %>% 
  mutate(metastatisi = ifelse(grepl('YES', metastatic_disease_confirmed), 1, 0)) %>%
  mutate(gender = ifelse(grepl('FEMALE', gender), 1, 0)) %>%
  dplyr::select(1,3,6,23)

traits_ss <- traits_ss %>% mutate(status = ifelse(grepl('Alive', vital_status), 1, 0)) %>%
  mutate(otm = ifelse(grepl('Dead', history_other_malignancy), 1, 0))

traits_ss <- traits_ss %>% dplyr::select(1,4:6)


# binarize categorical variables


metastasis_site <- binarizeCategoricalColumns(col_metada_ss$metastasis_site,
                                              includePairwise = FALSE,
                                              includeLevelVsAll = TRUE,
                                              minCount = 1)

ethnicity <- binarizeCategoricalColumns(col_metada_ss$ethnicity,
                                   includePairwise = FALSE,
                                   includeLevelVsAll = TRUE,
                                   minCount = 1)

traits_ss <- cbind(traits_ss, metastasis_site, ethnicity)

# Define numbers of genes and samples
nSamples_ss <- nrow(norm.counts_ss)
nGenes_ss <- ncol(norm.counts_ss)


module.trait.corr_ss <- cor(module_eigengenes_ss, traits_ss, use = 'p')
module.trait.corr.pvals_ss <- corPvalueStudent(module.trait.corr_ss, nSamples_ss)



# visualize module-trait association as a heatmap

heatmap.data_ss <- merge(module_eigengenes_ss, traits_ss, by = 'row.names')

head(heatmap.data_ss)

heatmap.data_ss <- heatmap.data_ss[,-1]
colnames(heatmap.data_ss)[63:65] <- c("Gender","Metastasis","Status")
colnames(heatmap.data_ss)[67] <- c("Lung") 

pdf("ss_wgcna.pdf")
CorLevelPlot(heatmap.data_ss,
             x = names(heatmap.data_ss)[c(63:65,67)],
             y = names(heatmap.data_ss)[1:62],
             col = c("blue1", "skyblue", "white", "pink", "red"),
             cexLabY = 0.4, cexCorval = 0.4
             )

dev.off()

#Vizualize genes in a module
module.gene.mapping <- as.data.frame(bwnet_ss$colors)
module.gene.mapping %>% 
  filter(`bwnet_ss$colors` == 'coral1') %>% 
  rownames()

#Indentify Hub Genes
#Metastasis
ss_m <- as.data.frame(traits_ss$Metastasis)
names(ss_m) <- "ss_m"

#names (colors) of the modules
modNames_ss <-  substring(names(module_eigengenes_ss), 3)

geneModuleMembership_ss <- as.data.frame(cor(norm.counts_ss, module_eigengenes_ss, use = "p"))
MMPvalue_ss <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_ss), nSamples_ss))


names(geneModuleMembership_ss) <- paste("MM", modNames_ss, sep="")
names(MMPvalue_ss) <- paste("p.MM", modNames_ss, sep="")

geneTraitSignificance_ss <- as.data.frame(cor(norm.counts_ss, ss_m, use = "p"))
GSPvalue_ss <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_ss), nSamples_ss))
GSPvalue_ss_adjust <- as.data.frame(p.adjust(as.matrix(GSPvalue_ss), method = "fdr"))

names(geneTraitSignificance_ss) <- paste("GS.", names(ss_m), sep="")
names(GSPvalue_ss) <- paste("p.GS.", names(ss_m), sep="")
names(GSPvalue_ss_adjust) <- paste("p.Adjuts.GS.", names(ss_m), sep="")

module_ss <- "coral1" #########################putting the color below the plot
column_ss <- match(module_ss, modNames_ss)
moduleGenes_ss <- moduleColors_ss==module_ss

verboseScatterplot(abs(geneModuleMembership_ss[moduleGenes_ss, column_ss]),
                   abs(geneTraitSignificance_ss[moduleGenes_ss, 1]),
                   xlab = paste("Module Membership in", module_ss, "module"),
                   ylab = "Gene significance for ss",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_ss)

geneInfo0_ss <- data.frame(EST = colnames(norm.counts_ss),
                          moduleColor = moduleColors_ss,
                          geneTraitSignificance_ss,
                          GSPvalue_ss,
                          GSPvalue_ss_adjust)


modOrder_ss <- order(-abs(cor(module_eigengenes_ss, ss_m, use = "p")))
for (mod in 1:ncol(geneModuleMembership_ss))
{
  oldNames = names(geneInfo0_ss)
  geneInfo0_ss = data.frame(geneInfo0_ss, geneModuleMembership_ss[, modOrder_ss[mod]], 
                            MMPvalue_ss[, modOrder_ss[mod]]);
  names(geneInfo0_ss) = c(oldNames, paste("MM.", modNames_ss[modOrder_ss[mod]], sep=""),
                          paste("p.MM.", modNames_ss[modOrder_ss[mod]], sep=""))
}
geneOrder_ss <- order(geneInfo0_ss$moduleColor, -abs(geneInfo0_ss$GS.ss_m))
geneInfo_ss <- geneInfo0_ss[geneOrder_ss, ]

coral_ss <- geneInfo_ss %>%
  filter(moduleColor=="coral1") %>%
  filter(abs(GS.ss_m) >= 0.2) %>%
  filter(p.Adjuts.GS.ss_m <= 0.05) %>%
  filter(abs(MM.coral1) >= 0.8) %>%
  filter(p.MM.coral1<= 0.05)

#Status
ss_s <- as.data.frame(traits_ss$Status)
names(ss_s) <- "ss_s"

#names (colors) of the modules
modNames_ss <-  substring(names(module_eigengenes_ss), 3)

geneModuleMembership_ss <- as.data.frame(cor(norm.counts_ss, module_eigengenes_ss, use = "p"))
MMPvalue_ss <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_ss), nSamples_ss))


names(geneModuleMembership_ss) <- paste("MM", modNames_ss, sep="")
names(MMPvalue_ss) <- paste("p.MM", modNames_ss, sep="")

geneTraitSignificance_ss <- as.data.frame(cor(norm.counts_ss, ss_s, use = "p"))
GSPvalue_ss <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_ss), nSamples_ss))
GSPvalue_ss_adjust <- as.data.frame(p.adjust(as.matrix(GSPvalue_ss), method = "fdr"))

names(geneTraitSignificance_ss) <- paste("GS.", names(ss_s), sep="")
names(GSPvalue_ss) <- paste("p.GS.", names(ss_s), sep="")
names(GSPvalue_ss_adjust) <- paste("p.Adjuts.GS.", names(ss_s), sep="")

module_ss <- "lightsteelblue1" #########################putting the color below the plot
column_ss <- match(module_ss, modNames_ss)
moduleGenes_ss <- moduleColors_ss==module_ss

verboseScatterplot(abs(geneModuleMembership_ss[moduleGenes_ss, column_ss]),
                   abs(geneTraitSignificance_ss[moduleGenes_ss, 1]),
                   xlab = paste("Module Membership in", module_ss, "module"),
                   ylab = "Gene significance for ss",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_ss)

geneInfo0_ss <- data.frame(EST = colnames(norm.counts_ss),
                          moduleColor = moduleColors_ss,
                          geneTraitSignificance_ss,
                          GSPvalue_ss,
                          GSPvalue_ss_adjust)


modOrder_ss <- order(-abs(cor(module_eigengenes_ss, ss_s, use = "p")))
for (mod in 1:ncol(geneModuleMembership_ss))
{
  oldNames = names(geneInfo0_ss)
  geneInfo0_ss = data.frame(geneInfo0_ss, geneModuleMembership_ss[, modOrder_ss[mod]], 
                            MMPvalue_ss[, modOrder_ss[mod]]);
  names(geneInfo0_ss) = c(oldNames, paste("MM.", modNames_ss[modOrder_ss[mod]], sep=""),
                          paste("p.MM.", modNames_ss[modOrder_ss[mod]], sep=""))
}
geneOrder_ss <- order(geneInfo0_ss$moduleColor, -abs(geneInfo0_ss$GS.ss_s))
geneInfo_ss <- geneInfo0_ss[geneOrder_ss, ]

#Indetify Hub genes
ls_ss <- geneInfo_ss %>%
  filter(moduleColor=="lightsteelblue1") %>%
  filter(abs(GS.ss_s) >= 0.2) %>%
  filter(p.Adjuts.GS.ss_s <= 0.05) %>%
  filter(abs(MM.lightsteelblue1) >= 0.8) %>%
  filter(p.MM.lightsteelblue1<= 0.05)

#Gender
ss_g <- as.data.frame(traits_ss$Gender)
names(ss_g) <- "ss_g"

#names (colors) of the modules
modNames_ss <-  substring(names(module_eigengenes_gs), 3)

geneModuleMembership_ss <- as.data.frame(cor(norm.counts_ss, module_eigengenes_ss, use = "p"))
MMPvalue_ss <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_ss), nSamples_ss))


names(geneModuleMembership_ss) <- paste("MM", modNames_ss, sep="")
names(MMPvalue_ss) <- paste("p.MM", modNames_ss, sep="")

geneTraitSignificance_ss <- as.data.frame(cor(norm.counts_ss, ss_g, use = "p"))
GSPvalue_ss <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_ss), nSamples_ss))
GSPvalue_ss_adjust <- as.data.frame(p.adjust(as.matrix(GSPvalue_ss), method = "fdr"))

names(geneTraitSignificance_ss) <- paste("GS.", names(ss_g), sep="")
names(GSPvalue_ss) <- paste("p.GS.", names(ss_g), sep="")
names(GSPvalue_ss_adjust) <- paste("p.Adjuts.GS.", names(ss_g), sep="")

module_ss <- "darkgreen" #########################putting the color below the plot
column_ss <- match(module_ss, modNames_ss)
moduleGenes_ss <- moduleColors_ss==module_ss

verboseScatterplot(abs(geneModuleMembership_ss[moduleGenes_ss, column_ss]),
                   abs(geneTraitSignificance_ss[moduleGenes_ss, 1]),
                   xlab = paste("Module Membership in", module_ss, "module"),
                   ylab = "Gene significance for ss",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_ss)

geneInfo0_ss <- data.frame(EST = colnames(norm.counts_ss),
                          moduleColor = moduleColors_ss,
                          geneTraitSignificance_ss,
                          GSPvalue_ss,
                          GSPvalue_ss_adjust)


modOrder_ss <- order(-abs(cor(module_eigengenes_ss, ss_s, use = "p")))
for (mod in 1:ncol(geneModuleMembership_ss))
{
  oldNames = names(geneInfo0_ss)
  geneInfo0_ss = data.frame(geneInfo0_ss, geneModuleMembership_ss[, modOrder_ss[mod]], 
                            MMPvalue_ss[, modOrder_ss[mod]]);
  names(geneInfo0_ss) = c(oldNames, paste("MM.", modNames_ss[modOrder_ss[mod]], sep=""),
                          paste("p.MM.", modNames_ss[modOrder_ss[mod]], sep=""))
}
geneOrder_ss <- order(geneInfo0_ss$moduleColor, -abs(geneInfo0_ss$GS.ss_s))
geneInfo_ss <- geneInfo0_ss[geneOrder_ss, ]

#Indetify Hub genes
dg_ss <- geneInfo_ss %>%
  filter(moduleColor=="darkgreen") %>%
  filter(abs(GS.ss_s) >= 0.2) %>%
  filter(p.Adjuts.GS.ss_s <= 0.05) %>%
  filter(abs(MM.darkgreen) >= 0.8) %>%
  filter(p.MM.darkgreen <= 0.05)

####Save Hub genes
write.xlsx(coral_ss, "wgcna_ss_coral.xlsx")
write.xlsx(ls_ss, "wgcna_ss_ls.xlsx")
write.xlsx(dg_ss, "wgcna_ss_dg.xlsx")

