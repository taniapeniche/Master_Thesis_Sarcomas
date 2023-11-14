#CommonData----
#Survival Data
time<-as.data.frame(CI_TCGA_mRNA[,c("death_days_to","last_contact_days_to")])
time[time== "[Not Applicable]" ] <- NA
time[time== "[Not Available]" ] <- NA
time[time== "-Inf" ] <- NA
time[,1]<-ifelse(is.na(time[,1]), time[,2], time[,1])
time[,2]<-NULL
time$death_days_to<-as.numeric(time$death_days_to)/365
time[is.na(time)] <- 0

vitalStatus <- as.data.frame(CI_TCGA_mRNA[,"vital_status"])
vitalStatus[vitalStatus== "Alive"] <- 1
vitalStatus[vitalStatus== "Dead" ] <- 0

metastasis <- as.data.frame(CI_TCGA_mRNA[,"metastatic_disease_confirmed"])
metastasis[metastasis== "YES"] <- 1
metastasis[metastasis== "NO" ] <- 0

data <- cbind(time, vitalStatus, metastasis)
colnames(data)<-c("time", "status", "metastasis")
data$time<-as.numeric(data$time)
data$status <-as.numeric(data$status)
data$gender <- ifelse(CI_TCGA_mRNA$gender=="FEMALE", 1, 0)

#common----
common_genes <- blue$EST
both <- intersect(rownames(CI_TCGA_mRNA),names(mRNA_fpkm))
data_common <- data[both,]
surv_common<-Surv(data_common$time,data_common$status)

for(i in common_genes){
  geneexp_status <- t(as.data.frame(ifelse(mRNA_fpkm[i,] > median(as.numeric(mRNA_fpkm[i,])), "High", "Low")))
  data_common[,i] <- geneexp_status
}

surv_fit_common <- vector("list", length=length(common_genes))
names(surv_fit_common) <- common_genes
surv_dif_common <- vector("list", length=length(common_genes))
names(surv_dif_common) <- common_genes
pValue_cc<-c()

for(i in common_genes){
  surv_fit <- survfit(as.formula(paste("surv_common ~", i)), data = data_common)
  surv_fit_common[[i]] <- surv_fit
  surv_dif <- survdiff(as.formula(paste("surv_common ~", i)), data = data_common)
  surv_dif_common[[i]] <- surv_dif
  pValue_mp <- pchisq(surv_dif$chisq, length(surv_dif$n)-1, lower.tail = FALSE)
  pValue_cc<-c(pValue_cc,as.numeric(pValue_mp))
}

pValue_cc_adjust<-p.adjust(pValue_cc, method="fdr")

CC_table <- as.data.frame(cbind(pValue_cc,pValue_cc_adjust))
row.names(CC_table) <- common_genes

CC_table %>% filter(pValue_cc_adjust <= 0.05)

#Leiomyosarcoma----
both <- intersect(rownames(CI_TCGA_mRNA),names(ts_geneexp[["LMS"]]))
data_lms <- data[both,]
geneExpStatus_lms <- t(as.data.frame(ifelse(lms_mrna_fpkm["CD74",] > median(as.numeric(lms_mrna_fpkm["CD74",])), "High", "Low")))
data_lms <- cbind(data_lms,geneExpStatus_lms)

surv_lms<-Surv(data_lms$time,data_lms$status)
surv_fit_lms<- survfit(surv_lms ~ CD74, data=data_lms)

surv_groups_lms <- survdiff(surv_lms ~ CD74, data=data_lms)
pValue_sv_lms <- pchisq(surv_groups_lms$chisq, length(surv_groups_lms$n)-1, lower.tail = FALSE)
pValue_sv_lms<-as.numeric(pValue_sv_lms)

ggsurvplot(surv_fit_lms, conf.int = F,risk.table = TRUE,pval = TRUE, 
           p.adjust = TRUE, p.adjust.methods="FDR", xlab = "Time (years)", 
           title = "CD74_lms" ,censor = T, data=data_lms)


#MPNST----
mpnst_genes <- c(maroon_mpnst$EST, dg_mpnst$EST, salmon4_mpnst$EST,
                 mb_mpnst$EST, black_mpnst$EST)



both <- intersect(rownames(CI_TCGA_mRNA),names(ts_geneexp[["MPNST"]]))
data_mpnst <- data[both,]
surv_mpnst<-Surv(data_mpnst$time,data_mpnst$status)

for(i in mpnst_genes){
  geneexp_status <- t(as.data.frame(ifelse(mpnst_mrna_fpkm[i,] > median(as.numeric(mpnst_mrna_fpkm[i,])), "High", "Low")))
  data_mpnst[,i] <- geneexp_status
}

surv_fit_mpnst <- vector("list", length=length(mpnst_genes))
names(surv_fit_mpnst) <- mpnst_genes
surv_dif_mpnst <- vector("list", length=length(mpnst_genes))
names(surv_dif_mpnst) <- mpnst_genes
pValue_sv<-c()

for(i in mpnst_genes){
  surv_fit <- survfit(as.formula(paste("surv_mpnst ~", i)), data = data_mpnst)
  surv_fit_mpnst[[i]] <- surv_fit
  surv_dif <- survdiff(as.formula(paste("surv_mpnst ~", i)), data = data_mpnst)
  surv_dif_mpnst[[i]] <- surv_dif
  pValue_mp <- pchisq(surv_dif$chisq, length(surv_dif$n)-1, lower.tail = FALSE)
  pValue_sv<-c(pValue_sv,as.numeric(pValue_mp))
}

pValue_sv_adjust<-p.adjust(pValue_sv, method="fdr")

SV_table <- as.data.frame(cbind(pValue_sv,pValue_sv_adjust))
row.names(SV_table) <- mpnst_genes

SV_table %>% filter(pValue_sv_adjust <= 0.05)



#SS----
ss_genes <- c(coral_ss$EST,ls_ss$EST,dg_ss$EST)
both <- intersect(rownames(CI_TCGA_mRNA),names(ts_geneexp[["SS"]]))
data_ss <- data[both,]
surv_ss<-Surv(data_ss$time,data_ss$status)

for(i in ss_genes){
  geneexp_status <- t(as.data.frame(ifelse(ss_mrna_fpkm[i,] > median(as.numeric(ss_mrna_fpkm[i,])), "High", "Low")))
  data_ss[,i] <- geneexp_status
}

surv_fit_ss <- vector("list", length=length(ss_genes))
names(surv_fit_ss) <- ss_genes
surv_dif_ss <- vector("list", length=length(ss_genes))
names(surv_dif_ss) <- ss_genes
pValue_ss<-c()

for(i in ss_genes){
  surv_fit <- survfit(as.formula(paste("surv_ss ~", i)), data = data_ss)
  surv_fit_ss[[i]] <- surv_fit
  surv_dif <- survdiff(as.formula(paste("surv_ss ~", i)), data = data_ss)
  surv_dif_ss[[i]] <- surv_dif
  pValue_mp <- pchisq(surv_dif$chisq, length(surv_dif$n)-1, lower.tail = FALSE)
  pValue_ss<-c(pValue_ss,as.numeric(pValue_mp))
}

pValue_ss_adjust<-p.adjust(pValue_ss, method="fdr")

SS_table <- as.data.frame(cbind(pValue_ss,pValue_ss_adjust))
row.names(SS_table) <- ss_genes

SS_table %>% filter(pValue_ss_adjust <= 0.05)