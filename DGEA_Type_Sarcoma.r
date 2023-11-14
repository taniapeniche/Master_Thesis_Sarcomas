#Data_Preparation----
counts_data <- geneexptum
clinical_info <- CI_TCGA

#sample with both data
id_common <- intersect(colnames(counts_data), clinical_info$Tumor_Sample_Barcode)

clinical_info <- clinical_info %>% filter(Tumor_Sample_Barcode %in% id_common)
rownames(clinical_info) <- clinical_info$Tumor_Sample_Barcode

counts_data <- counts_data[,rownames(clinical_info)]

samples <- colnames(counts_data)
condition_ts <- c()

for(i in 1:length(samples)){
  for(j in sarcomatypes){
    if(samples[i] %in% names(ts_geneexp[[j]])){
      condition_ts <- c(condition_ts, j)
    }
  }
}

clinical_info$TS <- condition_ts


#####PCA####
pca.res <- pca(counts_data, metadata=clinical_info)

#Plot variance explained by each component
screeplot(pca.res)

#Plot 2 selected components/eigenvectors
biplot(pca.res)
biplot(pca.res, colby = "TS", hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
biplot(pca.res, lab="", colby="TS", hline = 0, vline = 0,legendPosition = 'top') # Biplot without sample names
biplot(pca.res, x="PC1", y="PC3",lab="",colby="TS", hline = 0, vline = 0,legendPosition = 'top') # Biplot with PC1 and PC3


#Plot several components
pairsplot(pca.res, colby="TS")


#DGE----
y_ts <- DGEList(counts=counts_data, group = condition_ts)

keep <- filterByExpr(y_ts, group=condition_ts)
y_ts <- y_ts[keep, keep.lib.sizes=FALSE]

# Normalize for library sizes
y_ts <- calcNormFactors(y_ts)
y_ts$samples

#####Design Matrix and Dispersion####
design_matrix <- model.matrix(~0+condition_ts)
colnames(design_matrix) <- gsub("condition_ts", "", colnames(design_matrix))
rownames(design_matrix) <- colnames(counts_data)

y_ts <- estimateDisp(y_ts, design_matrix)
plotBCV(y_ts)


fit <- glmQLFit(y_ts,design_matrix)

names_matrix <- colnames(design_matrix)

#Contrast Matrix and QLF test---
####DL###
contrast_matrix_1.1 <- makeContrasts(DL-LMS, levels = design_matrix)
qlf_1.1 <- glmQLFTest(fit, contrast =  contrast_matrix_1.1)
contrast_matrix_1.2 <- makeContrasts(DL-MPNST, levels = design_matrix)
qlf_1.2 <- glmQLFTest(fit, contrast =  contrast_matrix_1.2)
contrast_matrix_1.3 <- makeContrasts(DL-SS, levels = design_matrix)
qlf_1.3 <- glmQLFTest(fit, contrast =  contrast_matrix_1.3)
contrast_matrix_1.4 <- makeContrasts(DL-UPS, levels = design_matrix)
qlf_1.4 <- glmQLFTest(fit, contrast =  contrast_matrix_1.4)
contrast_matrix_1.5 <- makeContrasts(DL-MF, levels = design_matrix)
qlf_1.5 <- glmQLFTest(fit, contrast =  contrast_matrix_1.5)

####LMS###
contrast_matrix_2.1 <- makeContrasts(LMS-DL, levels = design_matrix)
qlf_2.1 <- glmQLFTest(fit, contrast =  contrast_matrix_2.1)
contrast_matrix_2.2 <- makeContrasts(LMS-MPNST, levels = design_matrix)
qlf_2.2<- glmQLFTest(fit, contrast =  contrast_matrix_2.2)
contrast_matrix_2.3 <- makeContrasts(LMS-SS, levels = design_matrix)
qlf_2.3 <- glmQLFTest(fit, contrast =  contrast_matrix_2.3)
contrast_matrix_2.4 <- makeContrasts(LMS-UPS, levels = design_matrix)
qlf_2.4 <- glmQLFTest(fit, contrast =  contrast_matrix_2.4)
contrast_matrix_2.5 <- makeContrasts(LMS-MF, levels = design_matrix)
qlf_2.5 <- glmQLFTest(fit, contrast =  contrast_matrix_2.5)

####MPNST####
contrast_matrix_3.1 <- makeContrasts(MPNST-DL, levels = design_matrix)
qlf_3.1 <- glmQLFTest(fit, contrast =  contrast_matrix_3.1)
contrast_matrix_3.2 <- makeContrasts(MPNST-LMS, levels = design_matrix)
qlf_3.2 <- glmQLFTest(fit, contrast =  contrast_matrix_3.2)
contrast_matrix_3.3 <- makeContrasts(MPNST-SS, levels = design_matrix)
qlf_3.3 <- glmQLFTest(fit, contrast =  contrast_matrix_3.3)
contrast_matrix_3.4 <- makeContrasts(MPNST-UPS, levels = design_matrix)
qlf_3.4 <- glmQLFTest(fit, contrast =  contrast_matrix_3.4)
contrast_matrix_3.5 <- makeContrasts(MPNST-MF, levels = design_matrix)
qlf_3.5 <- glmQLFTest(fit, contrast =  contrast_matrix_3.5)

#####SS####
contrast_matrix_4.1 <- makeContrasts(SS-DL, levels = design_matrix)
qlf_4.1 <- glmQLFTest(fit, contrast =  contrast_matrix_4.1)
contrast_matrix_4.2 <- makeContrasts(SS-LMS, levels = design_matrix)
qlf_4.2 <- glmQLFTest(fit, contrast =  contrast_matrix_4.2)
contrast_matrix_4.3 <- makeContrasts(SS-MPNST, levels = design_matrix)
qlf_4.3 <- glmQLFTest(fit, contrast =  contrast_matrix_4.3)
contrast_matrix_4.4 <- makeContrasts(SS-UPS, levels = design_matrix)
qlf_4.4 <- glmQLFTest(fit, contrast =  contrast_matrix_4.4)
contrast_matrix_4.5 <- makeContrasts(SS-MF, levels = design_matrix)
qlf_4.5 <- glmQLFTest(fit, contrast =  contrast_matrix_4.5)

####UPS###
contrast_matrix_5.1 <- makeContrasts(UPS-DL, levels = design_matrix)
qlf_5.1 <- glmQLFTest(fit, contrast = contrast_matrix_5.1)
contrast_matrix_5.2 <- makeContrasts(UPS-LMS, levels = design_matrix)
qlf_5.2 <- glmQLFTest(fit, contrast =  contrast_matrix_5.2)
contrast_matrix_5.3 <- makeContrasts(UPS-MPNST, levels = design_matrix)
qlf_5.3 <- glmQLFTest(fit, contrast =  contrast_matrix_5.3)
contrast_matrix_5.4 <- makeContrasts(UPS-SS, levels = design_matrix)
qlf_5.4 <- glmQLFTest(fit, contrast =  contrast_matrix_5.4)
contrast_matrix_5.5 <- makeContrasts(UPS-MF, levels = design_matrix)
qlf_5.5 <- glmQLFTest(fit, contrast =  contrast_matrix_5.5)

####MF####
contrast_matrix_6.1 <- makeContrasts(MF-DL, levels = design_matrix)
qlf_6.1 <- glmQLFTest(fit, contrast = contrast_matrix_6.1)
contrast_matrix_6.2 <- makeContrasts(MF-LMS, levels = design_matrix)
qlf_6.2 <- glmQLFTest(fit, contrast =  contrast_matrix_6.2)
contrast_matrix_6.3 <- makeContrasts(MF-MPNST, levels = design_matrix)
qlf_6.3 <- glmQLFTest(fit, contrast =  contrast_matrix_6.3)
contrast_matrix_6.4 <- makeContrasts(MF-SS, levels = design_matrix)
qlf_6.4 <- glmQLFTest(fit, contrast =  contrast_matrix_6.4)
contrast_matrix_6.5 <- makeContrasts(MF-UPS, levels = design_matrix)
qlf_6.5 <- glmQLFTest(fit, contrast =  contrast_matrix_6.5)

####qlf
qlf_all <- list(qlf_1.1 = qlf_1.1, 
                qlf_1.2 = qlf_1.2, 
                qlf_1.3 = qlf_1.3,
                qlf_1.4 = qlf_1.4, 
                qlf_1.5 = qlf_1.5,
                qlf_2.1 = qlf_2.1, 
                qlf_2.2 = qlf_2.2,
                qlf_2.3 = qlf_2.3,
                qlf_2.4 = qlf_2.4,
                qlf_2.5 = qlf_2.5,
                qlf_3.1 = qlf_3.1, 
                qlf_3.2 = qlf_3.2,
                qlf_3.3 = qlf_3.3,
                qlf_3.4 = qlf_3.4,
                qlf_3.5 = qlf_3.5,
                qlf_4.1 = qlf_4.1,
                qlf_4.2 = qlf_4.2,
                qlf_4.3 = qlf_4.3,
                qlf_4.4 = qlf_4.4,
                qlf_4.5 = qlf_4.5,
                qlf_5.1 = qlf_5.1,
                qlf_5.2 = qlf_5.2,
                qlf_5.3 = qlf_5.3, 
                qlf_5.4 = qlf_5.4,
                qlf_5.5 = qlf_5.5,
                qlf_6.1 = qlf_6.1,
                qlf_6.2 = qlf_6.2,
                qlf_6.3 = qlf_6.3, 
                qlf_6.4 = qlf_6.4,
                qlf_6.5 = qlf_6.5)

#####Plots####
degs_all <- vector("list", length = length(qlf_all))
upReg_all <- vector("list", length = length(qlf_all))
downReg_all <- vector("list", length = length(qlf_all))


for(i in 1:length(qlf_all)){
  qlf <- qlf_all[[i]]
  DEGs <- topTags(qlf, nrow(qlf$table), p.value=0.05, adjust.method = "fdr")$table
  degs_all[[i]] <- DEGs
  upReg <- rownames(DEGs)[which(DEGs$logFC > 1)]
  upReg_all[[i]] <- upReg
  downReg <- rownames(DEGs)[which(DEGs$logFC < -1)]
  downReg_all[[i]] <- downReg
  allGenes <- topTags(qlf, n = nrow(qlf$table), p.value = 1)$table
  plotData <- cbind(allGenes$logFC, -log10(allGenes$FDR)); rownames(plotData) <- rownames(allGenes)
  plot(plotData, pch=20,col="gray",xlab="Biological Variation (log2 Fold-Change)", ylab="Statistical Significance (-log10 P-Value)")
  abline(h=-log10(0.05), v=c(-1,1),lty=2)
  points(plotData[upReg,], col="#B40041", pch=20)
  points(plotData[downReg,], col="#0BCACA", pch=20)
  if(length(upReg) > 0){
    if (length(upReg) >= 5){
      text(plotData[upReg[1:5],], labels=upReg[1:5],col="#B40041", pos=sample(c(1:3), size=8, replace=T), cex=0.7, offset = 1)
      }
    if (length(upReg) < 5){
      text(plotData[upReg[1:length(upReg)],], labels=upReg[1:length(upReg)],col="#B40041", pos=sample(c(1:3), size=8, replace=T), cex=0.7, offset = 1)
    }
  }
  if(length(downReg) > 0){
    if (length(downReg) >= 5){
      text(plotData[downReg[1:5],], labels=downReg[1:5],col="#0BCACA", pos=sample(c(1:3), size=10, replace=T), cex=0.8, offset = 1)
      }
    if (length(downReg) < 5){
      text(plotData[downReg[1:length(downReg)],], labels=downReg[1:length(downReg)],col="#0BCACA", pos=sample(c(1:3), size=10, replace=T), cex=0.8, offset = 1)
      }
  }
}

dev.off()

#####Genes differential expressed in all the comparations in each type of sarcoma####

degs_DL_up <- list(upReg_all[[1]], upReg_all[[2]], upReg_all[[3]], upReg_all[[4]], upReg_all[[5]])
degs_DL_up <- Reduce(intersect, degs_DL_up)

degs_DL_down <- list(downReg_all[[1]],downReg_all[[2]], downReg_all[[3]], downReg_all[[4]], downReg_all[[5]])
degs_DL_down <- Reduce(intersect, degs_DL_down)

degs_LMS_up <- list(upReg_all[[6]], upReg_all[[7]], upReg_all[[8]], upReg_all[[9]], upReg_all[[10]])
degs_LMS_up <- Reduce(intersect, degs_LMS_up)

degs_LMS_down <- list(downReg_all[[6]],downReg_all[[7]], downReg_all[[8]], downReg_all[[9]], downReg_all[[10]])
degs_LMS_down <- Reduce(intersect, degs_LMS_down)

degs_MPNST_up <- list(upReg_all[[11]], upReg_all[[12]], upReg_all[[13]], upReg_all[[14]], upReg_all[[15]])
degs_MPNST_up <- Reduce(intersect, degs_MPNST_up)

degs_MPNST_down <- list(downReg_all[[11]], downReg_all[[12]], downReg_all[[13]], downReg_all[[14]], downReg_all[[15]])
degs_MPNST_down <- Reduce(intersect, degs_MPNST_down)

degs_SS_up <-list(upReg_all[[16]], upReg_all[[17]], upReg_all[[18]], upReg_all[[19]], upReg_all[[20]])
degs_SS_up <- Reduce(intersect, degs_SS_up)

degs_SS_down <- list(downReg_all[[16]], downReg_all[[17]], downReg_all[[18]], downReg_all[[19]], downReg_all[[20]])
degs_SS_down <- Reduce(intersect, degs_SS_down)

degs_UPS_up <- list(upReg_all[[21]], upReg_all[[22]], upReg_all[[24]], upReg_all[[25]], upReg_all[[23]])
degs_UPS_up <- Reduce(intersect, degs_UPS_up)

degs_UPS_down <- list(downReg_all[[21]], downReg_all[[22]], downReg_all[[23]], downReg_all[[24]], downReg_all[[25]])
degs_UPS_down <- Reduce(intersect, degs_UPS_down)

degs_MF_up <- list(upReg_all[[26]], upReg_all[[27]], upReg_all[[28]], upReg_all[[29]], upReg_all[[30]])
degs_MF_up <- Reduce(intersect, degs_MF_up)

degs_MF_down <-  list(downReg_all[[26]], downReg_all[[27]], downReg_all[[28]], downReg_all[[29]], downReg_all[[30]])
degs_MF_down <- Reduce(intersect, degs_MF_down)

allgenes <- rownames(allGenes)

write(allgenes, "allgenes.txt")
write(degs_DL_up, "degs_DL_up.txt")
write(degs_DL_down, "degs_DL_down.txt")
write(degs_LMS_up, "degs_LMS_up.txt")
write(degs_LMS_down, "degs_LMS_down.txt")
write(degs_MPNST_up, "degs_MPNST_up.txt")
write(degs_MPNST_down, "degs_MPNST_down.txt")
write(degs_SS_up, "degs_SS_up.txt")
write(degs_SS_down, "degs_SS_down.txt")
write(degs_UPS_up, "degs_UPS_up.txt")
write(degs_UPS_down, "degs_UPS_down.txt")
write(degs_MF_up, "degs_MF_up.txt")
write(degs_MF_down, "degs_MF_down.txt")


#####plot individual####

#######Liposarcoma####
pdf("DGE_DL.pdf", height = 8, width = 8)

par(mfrow = c(3,2))

for(i in 1:5){
  qlf <- qlf_all[[i]]
  DEGs <- degs_all[[i]]
  upReg <- upReg_all[[i]] 
  downReg <- downReg_all[[i]]
  allGenes <- topTags(qlf, n = nrow(qlf$table), p.value = 1)$table
  plotData <- cbind(allGenes$logFC, -log10(allGenes$FDR)); rownames(plotData) <- rownames(allGenes)
  plot(plotData, pch=20,col="gray",xlab="Biological Variation (log2 Fold-Change)", ylab="Statistical Significance (-log10 P-Value)")
  abline(h=-log10(0.05), v=c(-1,1),lty=2)
  points(plotData[upReg,], col="#D8BFD8", pch=20)
  points(plotData[degs_DL_up,], col="#B40041", pch=20)
  points(plotData[downReg,], col="#AFEEEE", pch=20)
  points(plotData[degs_DL_down,], col="#000080", pch=20)
}

dev.off()

######Leiomyosarcoma####

pdf("DGE_LMS.pdf", height = 8, width = 8)

par(mfrow = c(3,2))

for(i in 6:10){
  qlf <- qlf_all[[i]]
  DEGs <- degs_all[[i]]
  upReg <- upReg_all[[i]] 
  downReg <- downReg_all[[i]]
  allGenes <- topTags(qlf, n = nrow(qlf$table), p.value = 1)$table
  plotData <- cbind(allGenes$logFC, -log10(allGenes$FDR)); rownames(plotData) <- rownames(allGenes)
  plot(plotData, pch=20,col="gray",xlab="Biological Variation (log2 Fold-Change)", ylab="Statistical Significance (-log10 P-Value)")
  abline(h=-log10(0.05), v=c(-1,1),lty=2)
  points(plotData[upReg,], col="#D8BFD8", pch=20)
  points(plotData[degs_LMS_up,], col="#B40041", pch=20)
  points(plotData[downReg,], col="#AFEEEE", pch=20)
  points(plotData[degs_LMS_down,], col="#000080", pch=20)
}

dev.off()

#####MPNST####

pdf("DGE_MPNST.pdf", height = 8, width = 8)

par(mfrow = c(3,2))

for(i in 11:15){
  qlf <- qlf_all[[i]]
  DEGs <- degs_all[[i]]
  upReg <- upReg_all[[i]] 
  downReg <- downReg_all[[i]]
  allGenes <- topTags(qlf, n = nrow(qlf$table), p.value = 1)$table
  plotData <- cbind(allGenes$logFC, -log10(allGenes$FDR)); rownames(plotData) <- rownames(allGenes)
  plot(plotData, pch=20,col="gray",xlab="Biological Variation (log2 Fold-Change)", ylab="Statistical Significance (-log10 P-Value)")
  abline(h=-log10(0.05), v=c(-1,1),lty=2)
  points(plotData[upReg,], col="#D8BFD8", pch=20)
  points(plotData[degs_MPNST_up,], col="#B40041", pch=20)
  points(plotData[downReg,], col="#AFEEEE", pch=20)
  points(plotData[degs_MPNST_down,], col="#000080", pch=20)
}

dev.off()

#####Sinovial Sarcoma####

pdf("DGE_SS.pdf", height = 8, width = 8)

par(mfrow = c(3,2))

for(i in 16:20){
  qlf <- qlf_all[[i]]
  DEGs <- degs_all[[i]]
  upReg <- upReg_all[[i]] 
  downReg <- downReg_all[[i]]
  allGenes <- topTags(qlf, n = nrow(qlf$table), p.value = 1)$table
  plotData <- cbind(allGenes$logFC, -log10(allGenes$FDR)); rownames(plotData) <- rownames(allGenes)
  plot(plotData, pch=20,col="gray",xlab="Biological Variation (log2 Fold-Change)", ylab="Statistical Significance (-log10 P-Value)")
  abline(h=-log10(0.05), v=c(-1,1),lty=2)
  points(plotData[upReg,], col="#D8BFD8", pch=20)
  points(plotData[degs_SS_up,], col="#B40041", pch=20)
  points(plotData[downReg,], col="#AFEEEE", pch=20)
  points(plotData[degs_SS_down,], col="#000080", pch=20)
}

dev.off()

#####UPS####

pdf("DGE_UPS.pdf", height = 8, width = 8)

par(mfrow = c(3,2))

for(i in 21:25){
  qlf <- qlf_all[[i]]
  DEGs <- degs_all[[i]]
  upReg <- upReg_all[[i]] 
  downReg <- downReg_all[[i]]
  allGenes <- topTags(qlf, n = nrow(qlf$table), p.value = 1)$table
  plotData <- cbind(allGenes$logFC, -log10(allGenes$FDR)); rownames(plotData) <- rownames(allGenes)
  plot(plotData, pch=20,col="gray",xlab="Biological Variation (log2 Fold-Change)", ylab="Statistical Significance (-log10 P-Value)")
  abline(h=-log10(0.05), v=c(-1,1),lty=2)
  points(plotData[upReg,], col="#D8BFD8", pch=20)
  points(plotData[degs_UPS_up,], col="#B40041", pch=20)
  points(plotData[downReg,], col="#AFEEEE", pch=20)
  points(plotData[degs_UPS_down,], col="#000080", pch=20)
}

dev.off()

#####Myxofibrossarcoma####

pdf("DGE_MF.pdf", height = 8, width = 8)

par(mfrow = c(3,2))

for(i in 26:30){
  qlf <- qlf_all[[i]]
  DEGs <- degs_all[[i]]
  upReg <- upReg_all[[i]] 
  downReg <- downReg_all[[i]]
  allGenes <- topTags(qlf, n = nrow(qlf$table), p.value = 1)$table
  plotData <- cbind(allGenes$logFC, -log10(allGenes$FDR)); rownames(plotData) <- rownames(allGenes)
  plot(plotData, pch=20,col="gray",xlab="Biological Variation (log2 Fold-Change)", ylab="Statistical Significance (-log10 P-Value)")
  abline(h=-log10(0.05), v=c(-1,1),lty=2)
  points(plotData[upReg,], col="#D8BFD8", pch=20)
  points(plotData[degs_MF_up,], col="#B40041", pch=20)
  points(plotData[downReg,], col="#AFEEEE", pch=20)
  points(plotData[degs_MF_down,], col="#000080", pch=20)
}

dev.off()




