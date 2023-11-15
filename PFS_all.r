#Prepare Data for PFS----
#DataFrame with relevant clinical dat
PSD <- as.data.frame(cbind(dtt_20$PATIENT_ID,as.numeric(dtt_20$STOP_DATE - dtt_20$START_DATE),
                           dtt_20$Treatment_TYPE, dtt_20$SUBTYPE, dtt_20$DRUG))
names(PSD) <- c("PATIENT_ID", "TIME","Treatment_Type", "Subtype", "Drug")
View(PSD)

#Types of therapies/Create anaylisis groups
VCT <- c("VEGF INHIBITOR","CHEMOTHERAPY","TYROSINE KINASE INHIBITOR")
VT <- c("VEGF INHIBITOR","TYROSINE KINASE INHIBITOR")
Chem <- "CHEMOTHERAPY"
VPT <- c("VEGF INHIBITOR","PDGFR-ALPHA BLOCKER","TYROSINE KINASE INHIBITOR")
APM <- "ANTI-PD-1 MAB"
VPTC <- c("VEGF INHIBITOR","PDGFR-ALPHA BLOCKER","TYROSINE KINASE INHIBITOR","CHEMOTHERAPY")

tot <- c("VCT","VT","Chem","VPT","APM","VPTC")

tt <- vector("list", length = length(tot))
names(tt) <- tot

tt[[1]] <- VCT
tt[[2]] <- VT
tt[[3]] <- Chem
tt[[4]] <- VPT
tt[[5]] <- APM
tt[[6]] <- VPTC

#PFS for all types of therapies on PLCG1----
PLCG1_mut <- mut_20[mut_20$Hugo_Symbol=="PLCG1","Tumor_Sample_Barcode"][1:27]
Mut_ID <- sapply(PLCG1_mut,function(x)unlist(strsplit(x,"_"))[3])

surv_groups_PFS<- vector("list", length=length(tot))
names(surv_groups_PFS) <- tot
pValue_all_PFS<-vector("list", length=length(tot))
names(pValue_all_PFS) <- tot
pValue_sv_PFS<-("list")
surv_groups_PFS<- vector("list", length=length(tot))
names(surv_groups_PFS) <- tot
pValue_all_PFS<-vector("list", length=length(tot))
names(pValue_all_PFS) <- tot
pValue_sv_PFS<-("list")

for(i in 1:length(tt)){
  #data frame with all the information
  t <- PSD[PSD$Subtype %in% tt[[i]],]
  ID <- sapply(t$PATIENT_ID,function(x)unlist(strsplit(x,"_"))[2])
  t$type <- ID
  #Indentify the mut samples
  for (j in 1:length(t$type)){
    if (t$type[j] %in% Mut_ID){t$type[j]<-"Mut"}
    else {t$type[j]<-"WT"}
  }
  data_PFS[[i]]<-t
  #Surv analyses
  surv_tt <- Surv(as.numeric(t$TIME))
  surv_PFS[[i]]<- surv_tt
  sfit_tt <- survfit(surv_tt ~ type, data=t)
  Surv_fit_tt <- sfit_tt
  gp<-ggsurvplot(sfit_tt, conf.int = F,risk.table = TRUE, xlab = "Days", censor = T,data=t, pval = T, title = tot[i])
  print(gp)
  surv_groups <- survdiff(surv_tt ~ type, data=data_PFS[[i]])
  surv_groups_PFS[[i]] <- surv_groups
  pValue_sv_PFS[i] <- pchisq(surv_groups$chisq, length(surv_groups$n)-1, lower.tail = FALSE)
  pValue_sv_PFS<-as.numeric(pValue_sv_PFS)
}

pValue_adj_PFS<-p.adjust(pValue_sv_PFS, method="fdr")
names_SURV_PFS<-names(data_PFS)

#Other_Genes----
#Dataframe complete
####Chemeotherapy####
Chem_PFS <- PSD[PSD$Subtype=="CHEMOTHERAPY",]
#ID_Patients
ID_Chem <- sapply(Chem_PFS$PATIENT_ID,function(x)unlist(strsplit(x,"_"))[2])
Chem_PFS$type <- ID_Chem

#####PFS_Chemeo####
#Genes mut
genes_mut <- c("KDR","PIK3CA","PTPRB", "LRP2","MUC16","PKHD1","TTN","RYR2","POT1")

ID_Mut <- vector("list", length = length(genes_mut))
names(ID_Mut) <- genes_mut

data_mut <- vector("list", length = length(genes_mut))
names(data_mut) <- genes_mut

data_per_gene<-vector("list", length = length(genes_mut))
names(data_per_gene) <- genes_mut

Surv_per_gene<-vector("list", length = length(genes_mut))
names(Surv_per_gene) <- genes_mut

Surv_fit_gene <- vector("list", length = length(genes_mut))
names(Surv_fit_gene) <- genes_mut

for(i in 1:length(ID_Mut)){
  #Select mut id for each gene and save
  id_mut <- mut_20[mut_20$Hugo_Symbol==genes_mut[i],"Tumor_Sample_Barcode"]
  id_mut_s <- sapply(id_mut,function(x)unlist(strsplit(x,"_"))[3])
  ID_Mut[[i]]<-id_mut_s
  #Create a data frame with only the mut 
  df <- Chem_PFS[Chem_PFS$type %in% id_mut_s,]
  df$type <- genes_mut[i]
  data_mut[[i]]<-df
  Wt_mut <- Chem_PFS
  for (j in 1:length(Wt_mut$type)){
    if (Wt_mut$type[j] %in% id_mut_s){Wt_mut$type[j]<-"Mut"}
    else {Wt_mut$type[j]<-"WT"}
  }
  data_per_gene[[i]]<-Wt_mut
  surv_gene <- Surv(as.numeric(Wt_mut$TIME))
  Surv_per_gene[[i]]<- surv_gene
  sfit_gene <- survfit(surv_gene ~ type, data=Wt_mut)
  Surv_fit_gene[[i]] <- sfit_gene
  gp<-ggsurvplot(sfit_gene, conf.int = F,risk.table = TRUE, xlab = "Days", censor = T,data=Wt_mut, pval = T, title = genes_mut[i])
  print(gp)
}

#####PFS_Angio_VPT####
#Genes mut
#genes_mut <- c("KDR","PIK3CA","PTPRB", "LRP2","MUC16","PKHD1","TTN","RYR2","POT1")

VPT_PFS <- PSD[PSD$Subtype %in% VPT,]
#ID_Patients
VPT_ID <- sapply(VPT_PFS$PATIENT_ID,function(x)unlist(strsplit(x,"_"))[2])
VPT_PFS$type <- VPT_ID

ID_Mut_V <- vector("list", length = length(genes_mut))
names(ID_Mut_V) <- genes_mut
data_mut_V <- vector("list", length = length(genes_mut))
names(data_mut_V) <- genes_mut
data_per_gene_V<-vector("list", length = length(genes_mut))
names(data_per_gene_V) <- genes_mut
Surv_per_gene_V<-vector("list", length = length(genes_mut))
names(Surv_per_gene_V) <- genes_mut
Surv_fit_gene_V <- vector("list", length = length(genes_mut))
names(Surv_fit_gene_V) <- genes_mut
surv_PFS_genes_V<- vector("list", length=length(genes_mut))
names(surv_PFS_genes_V) <- genes_mut
pValue_genes_PFS_V<-vector("list", length=length(genes_mut))
names(pValue_genes_PFS_V) <- genes_mut
pValue_sv_genes_V<-("list")

for(i in 1:length(ID_Mut_V)){
  #Select mut id for each gene and save
  id_mut <- mut_20[mut$Hugo_Symbol==genes_mut[i],"Tumor_Sample_Barcode"]
  id_mut_s <- sapply(id_mut,function(x)unlist(strsplit(x,"_"))[3])
  ID_Mut_V[[i]]<-id_mut_s
  #Create a data frame with only the mut 
  df <- VPT_PFS[VPT_PFS$type %in% id_mut_s,]
  df$type <- genes_mut[i]
  data_mut_V[[i]]<-df
  Wt_mut <- VPT_PFS
  for (j in 1:length(Wt_mut$type)){
    if (Wt_mut$type[j] %in% id_mut_s){Wt_mut$type[j]<-"Mut"}
    else {Wt_mut$type[j]<-"WT"}
  }
  data_per_gene_V[[i]]<-Wt_mut
  surv_gene <- Surv(as.numeric(Wt_mut$TIME))
  Surv_per_gene_V[[i]]<- surv_gene
  sfit_gene <- survfit(surv_gene ~ type, data=Wt_mut)
  Surv_fit_gene_V[[i]] <- sfit_gene
  gp<-ggsurvplot(sfit_gene, conf.int = F,risk.table = TRUE, xlab = "Days", censor = T,data=Wt_mut, pval = T, title = genes_mut[i])
  print(gp)
  surv_groups <- survdiff(Surv_per_gene_V[[i]] ~ type, data=data_per_gene_V[[i]])
  surv_PFS_genes_V[[i]]<-surv_groups
  pValue_sv_genes_V[i] <- pchisq(surv_groups$chisq, length(surv_groups$n)-1, lower.tail = FALSE)
  pValue_sv_genes_V<-as.numeric(pValue_sv_genes_V)

}

pValue_adj_genes_V<-p.adjust(pValue_sv_genes_V, method="fdr")
names_SURV_genes_V<-names(genes_mut)

#Createa a table
SV_genes_V<-data.frame(round(pValue_sv_genes_V, digits = 4),round(pValue_adj_genes_V, digits = 4))
row.names(SV_genes_V)<-names(genes_mut)
colnames(SV_genes_V)<-c("p-Value", "Adj. p-Value")

Table_genes_V<-SV_genes_V
row.names(Table_genes_V)<-genes_mut