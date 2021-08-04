library(TwoSampleMR)
library(dplyr)
'%ni%' <- Negate('%in%')
library(MRInstruments)
library(MVMR)

Pathway_geno="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_PUFA_07_16/SNP/"
Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_PUFA_07_16/result/"

COVID_LIST=c("HGI_round_4_A2","HGI_round_4_B1","HGI_round_4_B2","HGI_round_4_C2","HGI_round_5_A2_eur","HGI_round_5_B1_eur","HGI_round_5_B2_eur","HGI_round_5_C2_eur")
#COVID_LIST=c("")

My_MR <- function(exp_dat,outcome_dat) {
  rm(dat)
  try(dat<- exp_dat %>% inner_join(outcome_dat, by= "SNP"), silent=TRUE)
  if (exists("dat")==TRUE && length(dat$SNP)>0) {
    dat$beta.exposure=dat$beta.exposure/dat$SD
    dat$se.exposure=dat$se.exposure/dat$SD
    dat$mr_keep=TRUE
    
    if (length(which(dat$mr_keep=='TRUE'))>2) {
      #Horizontal pleiotropy
      pleiotropy=mr_pleiotropy_test(dat)
      Heterogeneity=mr_heterogeneity(dat)
      res <- mr(dat, method_list=c("mr_ivw_mre","mr_ivw","mr_egger_regression", "mr_two_sample_ml", "mr_simple_median", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
      IVW_MRE = res[res$method=='Inverse variance weighted (multiplicative random effects)',]
      MR_Egger = res[res$method=='MR Egger',]
      IVW = res[res$method=='Inverse variance weighted',]
      ml = res[res$method=='Maximum likelihood',]
      W_Med = res[res$method=='Weighted median',]
      W_Mod = res[res$method=='Weighted mode',]
      final_res <- IVW_MRE
    } else if (length(which(dat$mr_keep=='TRUE'))==1) {
      res <- mr(dat, method_list=c("mr_wald_ratio"))
      Wald = res[res$method=='Wald ratio',]
      final_res <- Wald
    } else if (length(which(dat$mr_keep=='TRUE'))==2) {
      res <- mr(dat, method_list=c("mr_ivw_mre","mr_ivw", "mr_two_sample_ml"))
      IVW_MRE = res[res$method=='Inverse variance weighted (multiplicative random effects)',]
      IVW = res[res$method=='Inverse variance weighted',]
      ml = res[res$method=='Maximum likelihood',]
      final_res <- IVW_MRE
    }
    
    final_res$Test=n
    names(final_res)[names(final_res) == "method"] <- "Methods"
    names(final_res)[names(final_res) == "nsnp"] <- "nsnps"
    names(final_res)[names(final_res) == "b"] <- "Beta"
    names(final_res)[names(final_res) == "se"] <- "SE"
    names(final_res)[names(final_res) == "exposure"] <- "trait"
    
    if (exists("Heterogeneity")==TRUE){
      Het_IVW=Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]
      names(Het_IVW)[names(Het_IVW) == "Q_pval"] <- "Het_IVW_pval"
      final_res<- Het_IVW %>% right_join(final_res, by= "id.exposure")
      Het_MR_Egger=Heterogeneity[Heterogeneity$method=='MR Egger',]
      names(Het_MR_Egger)[names(Het_MR_Egger) == "Q_pval"] <- "Het_Egger_pval"
      final_res<- Het_MR_Egger %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$Het_IVW_pval=NA
      final_res$Het_Egger_pval=NA
    }
    
    if (exists("pleiotropy")==TRUE){
      names(pleiotropy)[names(pleiotropy) == "egger_intercept"] <- "Egger_intercept"
      names(pleiotropy)[names(pleiotropy) == "pval"] <- "pval_intercept"
      final_res<- pleiotropy %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$Egger_intercept=NA
      final_res$pval_intercept=NA
    }
    
    if (exists("W_Mod")==TRUE){
      names(W_Mod)[names(W_Mod) == "b"] <- "b_W_Mod"
      names(W_Mod)[names(W_Mod) == "se"] <- "se_W_Mod"
      names(W_Mod)[names(W_Mod) == "pval"] <- "pval_W_Mod"
      final_res<- W_Mod %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_W_Mod=NA
      final_res$se_W_Mod=NA
      final_res$pval_W_Mod=NA
    }
    
    if (exists("W_Med")==TRUE){
      names(W_Med)[names(W_Med) == "b"] <- "b_W_Med"
      names(W_Med)[names(W_Med) == "se"] <- "se_W_Med"
      names(W_Med)[names(W_Med) == "pval"] <- "pval_W_Med"
      final_res<- W_Med %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_W_Med=NA
      final_res$se_W_Med=NA
      final_res$pval_W_Med=NA
    }
    
    if (exists("ml")==TRUE){
      names(ml)[names(ml) == "b"] <- "b_ml"
      names(ml)[names(ml) == "se"] <- "se_ml"
      names(ml)[names(ml) == "pval"] <- "pval_ml"
      final_res<- ml %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_ml=NA
      final_res$se_ml=NA
      final_res$pval_ml=NA
    }
    
    if (exists("Wald")==TRUE){
      names(Wald)[names(Wald) == "b"] <- "b_Wald"
      names(Wald)[names(Wald) == "se"] <- "se_Wald"
      names(Wald)[names(Wald) == "pval"] <- "pval_Wald"
      final_res<- Wald %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_Wald=NA
      final_res$se_Wald=NA
      final_res$pval_Wald=NA
    }
    
    if (exists("MR_Egger")==TRUE){
      names(MR_Egger)[names(MR_Egger) == "b"] <- "b_Egger"
      names(MR_Egger)[names(MR_Egger) == "se"] <- "se_Egger"
      names(MR_Egger)[names(MR_Egger) == "pval"] <- "pval_Egger"
      final_res<- MR_Egger %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_Egger=NA
      final_res$se_Egger=NA
      final_res$pval_Egger=NA
    }
    
    if (exists("IVW")==TRUE){
      names(IVW)[names(IVW) == "b"] <- "b_IVW"
      names(IVW)[names(IVW) == "se"] <- "se_IVW"
      names(IVW)[names(IVW) == "pval"] <- "pval_IVW"
      final_res<- IVW %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_IVW=NA
      final_res$se_IVW=NA
      final_res$pval_IVW=NA
    }
    
    if (exists("IVW_MRE")==TRUE){
      names(IVW_MRE)[names(IVW_MRE) == "b"] <- "b_IVW_MRE"
      names(IVW_MRE)[names(IVW_MRE) == "se"] <- "se_IVW_MRE"
      names(IVW_MRE)[names(IVW_MRE) == "pval"] <- "pval_IVW_MRE"
      final_res<- IVW_MRE %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_IVW_MRE=NA
      final_res$se_IVW_MRE=NA
      final_res$pval_IVW_MRE=NA
    }
    
    names(final_res)[names(final_res) == "id.exposure"] <- "id"
    
    final_res<-final_res%>%select(Test, id, trait, b_IVW_MRE, se_IVW_MRE, pval_IVW_MRE, b_Egger, se_Egger, pval_Egger,b_Wald, se_Wald, pval_Wald,
                                  Egger_intercept, pval_intercept, Het_IVW_pval, Het_Egger_pval,
                                  b_W_Med, se_W_Med, pval_W_Med, b_W_Mod, se_W_Mod, pval_W_Mod)
    
    final_res[is.na(final_res) ] <- "-"
    
    Outputfile=paste(Pathway_out,n, ".txt", sep="")
    write.table(final_res, file= Outputfile, col.names = F, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  }
}

for (n in COVID_LIST) {
  self_PLASMA_GWAS_id=c("ALA","LA", "GLA", "DGLA", "AA", "DPA-n3", "DHA")
  
  for (RBC_name in self_PLASMA_GWAS_id) {
    exp_dat_file=paste(Pathway_geno,"Plasma_",RBC_name,"_",n, ".txt", sep="")
    exp_dat1 <- read.csv(exp_dat_file,header=T, as.is=T,sep = "\t")
    exp_dat <- read_exposure_data(
      filename = exp_dat_file,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "beta.exposure",
      se_col = "se.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure",
      eaf_col = "eaf.exposure",
      pval_col = "pval.exposure",
      #id_col = "exposure",
      phenotype_col = "exposure",
      #units_col = "Units",
      gene_col = "Gene",
      samplesize_col = "samplesize.exposure"
    )
    
    exp_dat$id.exposure=paste("Plasma",  sep="")
    exp_dat$SD=exp_dat1$SD
    
    exp_dat <-clump_data(
      exp_dat,
      clump_kb = 10000,
      clump_r2 = 0.001,
      clump_p1 = 1,
      clump_p2 = 1,
      pop = "EUR"
    )
    
    outcomefile=paste(Pathway_geno,"Plasma_",RBC_name,"_",n, ".txt", sep="")
    ####### Change csv file
    outcome_dat=read_outcome_data(
      filename = outcomefile,
      snps = NULL,
      sep = "\t",
      #phenotype_col = "outcome",
      snp_col = "SNP",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      eaf_col = "eaf.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure",
      #samplesize_col = "all_meta_sample_N",
      pval_col = "pval.outcome")
    
    outcome_dat$outcome=n
    outcome_dat$id.outcome=n
    
    if (length(exp_dat$SNP)>0) {
      exp_dat=exp_dat[exp_dat$beta.exposure != 0, ]
      exp_dat=exp_dat[exp_dat$se.exposure != 0, ]
      if (length(exp_dat$SNP)>0) {
        My_MR(exp_dat,outcome_dat)
      }
    }
  }
  
  RBC_GWAS_id=c("ALA","LA", "GLA", "DGLA", "AA", "DPA-n3", "DTA")
  
  for (RBC_name in RBC_GWAS_id) {
    for (LD_num in c(0.001)) {
      exp_dat_file=paste(Pathway_geno,"RBC_",RBC_name,"_",n, ".txt", sep="")
      exp_dat1 <- read.csv(exp_dat_file,header=T, as.is=T,sep = "\t")
      exp_dat <- read_exposure_data(
        filename = exp_dat_file,
        sep = "\t",
        snp_col = "SNP",
        beta_col = "beta.exposure",
        se_col = "se.exposure",
        effect_allele_col = "effect_allele.exposure",
        other_allele_col = "other_allele.exposure",
        eaf_col = "eaf.exposure",
        pval_col = "pval.exposure",
        #id_col = "exposure",
        phenotype_col = "exposure",
        #units_col = "Units",
        gene_col = "Gene",
        samplesize_col = "samplesize.exposure"
      )
      
      exp_dat$id.exposure=paste("RBC",  sep="")
      exp_dat$SD=exp_dat1$SD
      
      exp_dat <-clump_data(
        exp_dat,
        clump_kb = 10000,
        clump_r2 = 0.001,
        clump_p1 = 1,
        clump_p2 = 1,
        pop = "EUR"
      )
      
      outcomefile=paste(Pathway_geno,"RBC_",RBC_name,"_",n, ".txt", sep="")
      ####### Change csv file
      outcome_dat=read_outcome_data(
        filename = outcomefile,
        snps = NULL,
        sep = "\t",
        #phenotype_col = "outcome",
        snp_col = "SNP",
        beta_col = "beta.outcome",
        se_col = "se.outcome",
        eaf_col = "eaf.exposure",
        effect_allele_col = "effect_allele.exposure",
        other_allele_col = "other_allele.exposure",
        #samplesize_col = "all_meta_sample_N",
        pval_col = "pval.outcome")
      
      outcome_dat$outcome=n
      outcome_dat$id.outcome=n
      
      if (length(exp_dat$SNP)>0) {
        exp_dat=exp_dat[exp_dat$beta.exposure != 0, ]
        exp_dat=exp_dat[exp_dat$se.exposure != 0, ]
        if (length(exp_dat$SNP)>0) {
          My_MR(exp_dat,outcome_dat)
        }
      }
    }
  }
}

