#####Create a list w/ the ID of the samples, a dataframe#####
Mut_genes <- unique(mut_20_b$Hugo_Symbol)
Mut_ID_all <- vector("list", length = length(Mut_genes))
names(Mut_ID_all) <- Mut_genes
Status_Cla <- vector("list", length = length(Mut_genes))
names(Status_Cla) <- Mut_genes
df_Status_all <- vector("list", length = length(Mut_genes))
names(df_Status_all) <- Mut_genes
WC_data <- vector("list", length = length(Mut_genes))
names(WC_data) <- Mut_genes
res_all<- vector("list", length = length(Mut_genes))
names(res_all) <- Mut_genes
genes_sign<-c()
p.Val_sign_genes <- c()
genes_sign_all<-c()
p.Val_sign_genes_all <- c()

for(i in 1:length(Mut_genes)){
  #Create Mut ID
  m <- mut_20_b[mut_20_b$Hugo_Symbol==Mut_genes[i],"Tumor_Sample_Barcode"][1:27]
  m_ID <- sapply(m,function(x)unlist(strsplit(x,"_"))[3])
  m_ID <- m_ID[!is.na(m_ID)]
  Mut_ID_all[[i]]<-m_ID
  #Classify as WT vs Mut
  TMB_Q <- ci_asp_20_b %>% dplyr::select(PATIENT_ID, Tumor_Sample_Barcode, TMB_NONSYNONYMOUS)
  a<-ifelse((TMB_Q$Tumor_Sample_Barcode %in% m), "Mut", "WT")
  t <- ifelse((TMB_Q$TMB_NONSYNONYMOUS >= 10), "High", "Low")
  Status_Cla[[i]]<-a
  #save the info info into the table
  TMB_Q$State <- a
  TMB_Q$TMB_G <- t
  TMB_Q$TMB_G <- as.factor(TMB_Q$TMB_G)
  df_Status_all[[i]]<-TMB_Q
  #Select the col that we want
  data<- TMB_Q  %>% dplyr::select(State, TMB_NONSYNONYMOUS)
  rownames(data)<-TMB_Q$Tumor_Sample_Barcode
  WC_data[[i]] <- data
  #Wilcox_test
  res <- wilcox.test(TMB_NONSYNONYMOUS ~ State , data = data, paired = F,conf.level = 0.99)
  res_all[[i]]<-res
  if(res$p.value<=0.0005){genes_sign<-c(genes_sign,Mut_genes[i])} 
  if(res$p.value<=0.0005){p.Val_sign_genes <- c(p.Val_sign_genes,round(res$p.value,7))}
  if(res$p.value<=0.05){genes_sign_all<-c(genes_sign_all,Mut_genes[i])} 
  if(res$p.value<=0.05){p.Val_sign_genes_all <- c(p.Val_sign_genes_all,round(res$p.value,7))}
}

genes_TMB_sign <- as.data.frame(cbind(genes_sign,p.Val_sign_genes))
genes_TMB_sign_all <- as.data.frame(cbind(genes_sign_all,p.Val_sign_genes_all))

#####Plot the significant ones####
Stat_all<- vector("list", length = length(genes_sign))
names(Stat_all) <- genes_sign
Cor_plot_all <- vector("list", length = length(genes_sign))
names(Cor_plot_all) <- genes_sign

for(i in genes_sign){
  data <- WC_data[[i]]
  stat.test <- data %>% 
    rstatix::wilcox_test(TMB_NONSYNONYMOUS ~ State) %>%
    add_significance()
  Stat_all[[i]]<- stat.test
}
  
#Plot
for(i in 1:length(Cor_plot_all)){
  data <- WC_data[[i]]
  stat.test <- Stat_all[[i]]
  stat.test <- stat.test %>% add_xy_position(x = "State")
  pl <- ggboxplot(data, x = "State", y = "TMB_NONSYNONYMOUS", 
                  color = "State", palette = c("#00AFBB", "#E7B800"),
                  ylab = "TMB", xlab = genes_sign[i]) + 
    stat_pvalue_manual(stat.test, tip.length = 0) +
    labs(subtitle = get_test_label(stat.test, detailed = TRUE))
  Cor_plot_all[[i]] <- pl
}

#Distrubution Analysis-----
#####Functions####
compare <- function(x,y){
  if(x>y){return(y)}
  else{return(x)}
}

selection <- function(x,y,z){
  if(length(x)>0 && length(y)>0 && length(z)>0){
    if(is.na(x)==F && is.na(y)==F && is.na(z)==F){
      a<-compare(x,y)
      b<-compare(a,z)
      return(b)
    }
    if(is.na(x)==F && is.na(y)==F && is.na(z)==T){
      a<-compare(x,y)
      return(a) 
    }
    if(is.na(x)==F && is.na(y)==T && is.na(z)==F){
      a<-compare(x,z)
      return(a)
    }
    if(is.na(x)==T && is.na(y)==F && is.na(z)==F){
      a<-compare(y,z)
      return(a)
    }
    if(is.na(x)==F && is.na(y)==T && is.na(z)==T){
      return(x)
    }
    if(is.na(x)==T && is.na(y)==F && is.na(z)==T){
      return(y)
    }
    if(is.na(x)==T && is.na(y)==T && is.na(z)==F){
      return(z)
    }
    if(is.na(x)==T && is.na(y)==T && is.na(z)==T){
      return(NA)
    }
  }
  else{
    if(length(x)==0 && length(y)>0 && length(z)>0){
      if(is.na(y)==T && is.na(z)==T){
        return("empty")
      }
      if(is.na(y)==F && is.na(z)==T){
        return(y)
      }
      if(is.na(y)==T && is.na(z)==F){
        return(z)
      }
      if(is.na(y)==F && is.na(z)==T){
        a<-compare(y,z)
        return(a)
      }
    }
    if(length(x)>0 && length(y)==0 && length(z)>0){
      if(is.na(x)==T && is.na(z)==T){
        return("empty")
      }
      if(is.na(x)==F && is.na(z)==T){
        return(x)
      }
      if(is.na(x)==T && is.na(z)==F){
        return(z)
      }
      if(is.na(x)==F && is.na(z)==T){
        a<-compare(x,z)
        return(a)
      }
    }
    if(length(x)>0 && length(y)>0 && length(z)==0){
      if(is.na(y)==T && is.na(x)==T){
        return("empty")
      }
      if(is.na(y)==F && is.na(x)==T){
        return(y)
      }
      if(is.na(y)==T && is.na(x)==F){
        return(x)
      }
      if(is.na(y)==F && is.na(x)==T){
        a<-compare(y,x)
        return(a)
      }
    }
    if(length(x)==0 && length(y)==0 && length(z)>0){
      return(z)
    }
    if(length(x)>0 && length(y)==0 && length(z)==0){
      return(x)
    }
    if(length(x)==0 && length(y)>0 && length(z)==0){
      return(y)
    }
    if(length(x)==0 && length(y)==0 && length(z)==0){
      return("empty")
    }
  }
}

#####Save Data####
load("~/Desktop/Sarcomas/Rfiles/TMB/TMB_Dist_Data.RData")

#####Type of distribution####
AIC_BIC <- vector("list",length = length(WC_data))
names(AIC_BIC) <-  names(WC_data)
stat <- c("WT","Mut")
gene <- names(AIC_BIC)
WT_Sel <- c()
MUT_Sel <- c()

for(i in 1:length(AIC_BIC)){
  data <- WC_data[[i]]
  WT <- c()
  MUT <- c()
  message_error <- c()
  for(j in stat){
    norm <- c("norm")
    lognorm <- c("lognorm")
    gamma <- c("gamma")
    data_s <- data %>% filter(State==j) %>% dplyr::select(TMB_NONSYNONYMOUS)
    AIC_norm <- c()
    BIC_norm <- c()
    AIC_lognorm <- c()
    BIC_lognorm <-c()
    AIC_gamma <- c()
    BIC_gamma <- c()
    tryCatch({
      fit_norm <- MASS::fitdistr(data_s$TMB_NONSYNONYMOUS, "normal")
      AIC_norm <- AIC(fit_norm)
      BIC_norm <- BIC(fit_norm)
      norm <- c(norm,AIC_norm, BIC_norm)
    }, error = function(e) {
        message_error <- c(message_error, "Error occurred: ", conditionMessage(e),names(AIC_BIC)[i])
    })
    tryCatch({
      fit_lognorm <- MASS::fitdistr(data_s$TMB_NONSYNONYMOUS, "lognormal")
      AIC_lognorm <- AIC(fit_lognorm)
      BIC_lognorm <- BIC(fit_lognorm)
      lognorm <- c(lognorm,AIC_lognorm, BIC_lognorm)
    }, error = function(e) {
      message_error <- c(message_error, "Error occurred: ", conditionMessage(e),names(AIC_BIC)[i])
    })
    tryCatch({
      fit_gamma<- MASS::fitdistr(data_s$TMB_NONSYNONYMOUS, "gamma")
      AIC_gamma<- AIC(fit_gamma)
      BIC_gamma<- BIC(fit_gamma)
      gamma <- c(gamma,AIC_gamma,BIC_gamma)
    }, error = function(e) {
      message_error <- c(message_error, "Error occurred: ", conditionMessage(e),names(AIC_BIC)[i])
    })
    if(j == "WT"){
      WT <-c(norm,lognorm,gamma,"WT") 
      val<-selection(AIC_norm,AIC_lognorm,AIC_gamma)
      if(length(val)>0){
        if(!(val=="empty")|is.na(val)==T){
          if(val==AIC_norm){WT_Sel<-c(WT_Sel,"norm")}
          if(val==AIC_lognorm){WT_Sel<-c(WT_Sel,"lognorm")}
          if(val==AIC_gamma){WT_Sel<-c(WT_Sel,"gamma")}
        }
        else{
          WT_Sel<-c(WT_Sel,val)
        }
      }
      else{
        WT_Sel<-c(WT_Sel,val)
      }
    }
    else{
      val<-selection(AIC_norm,AIC_lognorm,AIC_gamma)
      if(length(val)>0){
        if(!(val=="empty")|is.na(val)==T){
          if(val==AIC_norm){MUT_Sel<-c(MUT_Sel,"norm")}
          if(val==AIC_lognorm){MUT_Sel<-c(MUT_Sel,"lognorm")}
          if(val==AIC_gamma){MUT_Sel<-c(MUT_Sel,"gamma")}
        }
        else{
          MUT_Sel<-c(MUT_Sel,val)
        }
      }
      else{
        MUT_Sel<-c(MUT_Sel,val)
      }
    }
  }
  AIC_BIC[[i]]<-c(WT,MUT)
}


Dist_sel <- as.data.frame(cbind(gene,WT_Sel,MUT_Sel))

######Select significant genes#####

Dist_sel %>% group_by(MUT_Sel) %>% count()

Comp <- Dist_sel %>% filter(MUT_Sel %in% c("gamma","norm","Lognorm"))

gene_sel <- Comp$gene
ks_dist <- c()
gene_t <- c()

for(i in gene_sel){
  data <- WC_data[[i]]
  data_s <-  data %>% filter(State=="WT") %>% dplyr::select(TMB_NONSYNONYMOUS)
  data_m <- data %>% filter(State=="Mut") %>% dplyr::select(TMB_NONSYNONYMOUS)
  if(nrow(data_m)>1){
    ks<-ks.test(data_s$TMB_NONSYNONYMOUS,data_m$TMB_NONSYNONYMOUS)
    ks_dist<-c(ks_dist,ks$p.value)
    gene_t<-c(gene_t,i)
  }
}

tt_dist_adjust <- p.adjust(tt_dist,method = "fdr")
ks_dist_adjust <-p.adjust(ks_dist,method = "fdr")

tabel_sign <- as.data.frame(cbind(gene_t,ks_dist,ks_dist_adjust))


table_sign_dist_ks <- tabel_sign %>% filter(ks_dist_adjust <= 0.05)



#####Distribution Plots####
Dist_Plot <- vector("list",length = length(WC_data))
names(Dist_Plot) <- names(WC_data)

for(i in 1:length(Dist_Plot)){
  data <- WC_data[[i]]
  p <- data %>%
    ggplot( aes(x=TMB_NONSYNONYMOUS, fill=State)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 40) +
    scale_fill_manual(values=c("#77003A", "#E6CBFC")) +
    geom_density(alpha = .2, fill="#FF6666") +
    ggtitle(names(Dist_Plot)[i]) +
    theme_ipsum(base_family = "Helvetica", axis_text_size = 6) +
    labs(fill="")
  Dist_Plot[[i]]<-p
}

#Plot PLCG1
pdf("ks_test_PLCG1.pdf")
print(Dist_Plot[["PLCG1"]])
dev.off()
