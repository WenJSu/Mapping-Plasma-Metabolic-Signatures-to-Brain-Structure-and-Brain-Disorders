library(TwoSampleMR)
library(data.table)
library('ieugwasr')
library(RadialMR)

library(future)
library(future.apply)
plan(multisession, workers=parallel::detectCores()-100)

meta <- fread("NMR_Preprocessed.csv", sep = ",", header = TRUE)
meta_list <- colnames(meta)


IDP_df <- fread("brain_all_mergeID.csv", sep = ",", header = TRUE)

refer_df <- read.csv("output_all_pro2_reFDR_sig.csv")
meta_list_use <- unique(refer_df$meta)
colnames(refer_df)[colnames(refer_df) == "UKBID.x"] <- "UKBID"
refer_df = merge(IDP_df, refer_df, by = c('UKBID'))

bfile='EUR'

for(i in seq(1,length(meta_list_use))){
  metaname=meta_list_use[i]
  print(metaname)
  print("**************************")
  output1 <- paste0(i,"/",length(meta_list_use))
  print(output1)
  print("**************************")
  
  threshold=5e-8
  threshold_folder="MR_5e8"
  
  exposure_dat_0 <- read_exposure_data(
    filename = paste0(metaname,'/GWAS_merge.glm.linear'),
    sep = "\t",
    snp_col = "ID",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "OMITTED",
    eaf_col = "A1_FREQ",
    pval_col = "P"
    #log_pval = TRUE
  )
  exposure_dat_0$exposure <- rep(metaname,nrow(exposure_dat_0))
  # https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  exposure_dat_0$rsid <- exposure_dat_0$SNP
  exposure_dat_0$trait_id <- rep(metaname,nrow(exposure_dat_0))
  exposure_dat_0$pval <- exposure_dat_0$pval.exposure
  exposure_clumped <- ld_clump(
    exposure_dat_0,
    plink_bin = 'plink',
    bfile = bfile,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p = threshold,
    pop = "EUR"
  )
  
  IDP_df = refer_df[refer_df$meta==metaname,]
  
  for(j in seq(1,nrow(IDP_df))){

    output2 <- paste0(j,"/",nrow(IDP_df))
    print(output2)
    out = IDP_df$UKBID.y[j]
    IDP_num = IDP_df$Pheno[j]
    outputpath=paste0(threshold_folder,'/',metaname,'/',as.character(out),'/')
    filepath = paste0(outputpath, out, "_",metaname,"_sensitivity_MR_LOO")
    cnt = 1
    if (dir.exists(outputpath)) {
      message(paste("Directory already exists:", outputpath))
      if(file.exists(filepath)){
        file_content <- readLines(filepath)
        if(length(file_content)>0){
            next
        }else {
            cnt = 0 }

      }else{
        cnt = 0
      }
      
    } 
    
    if (cnt == 1) {
      dir.create(outputpath, recursive=TRUE)
      }
    
    print("......")
    meta_num <- paste0(i,"/",length(meta_list_use))
    print(meta_num)
    print("......")
    ## 1.Data sets preparation
    # Require two GWAS summary-level data with variants minor allele frequency (MAF) > 0.01; Removing the palindromic SNPs; Removing the long-range LD in the genome. (i.e. exposure_dat_0; outcome_dat_0)
    # can't done with plink, plink work for genotyping data, not summary data
    #########################################  2.MR Analysis preparation #########################################
    # a)Selection of conditional independent SNPs for exposures
    IDP_num_4 <- sprintf("%04d", IDP_num)
    outcomefile <- paste0(as.character(IDP_num_4),'.csv')
    
    snps=exposure_clumped$SNP
    outcome_dat_0 <- tryCatch({
      # 尝试执行第一段代码
      read_outcome_data(
        snps = snps,
        filename = outcomefile,
        sep = ",",
        snp_col = "rsid",
        beta_col = "beta",
        se_col = "se",
        effect_allele_col = "a2",
        other_allele_col = "a1",
        eaf_col = "af",  # merge列
        pval_col = "pval(-log10)",
        log_pval = TRUE
      )
    }, error = function(e) {
      # 如果第一段代码出错，执行备用代码
      message("第一段代码出错，正在执行备用代码...")
      read_outcome_data(
        snps = snps,
        filename = outcomefile,
        sep = ",",
        snp_col = "rsid",
        beta_col = "beta",
        se_col = "se",
        effect_allele_col = "a2",
        other_allele_col = "a1",
        eaf_col = "af",  # merge列
        pval_col = "pval..log10.",
        log_pval = TRUE
      )
    })
    
    
    outcome = as.character(out)
    outcome_dat_0$outcome <- rep(outcome,nrow(outcome_dat_0))
    sampleSize_outcome=243096
    outcome_dat_0$samplesize.outcome <- rep(sampleSize_outcome,nrow(outcome_dat_0))
    
    outcomefilename = outcome
    output1<-paste0(outputpath, outcomefilename,'_',metaname,'_QC_IVW') #output data IVW
    output2<-paste0(outputpath, outcomefilename,'_',metaname,'_QC_MR-egger') #output data MR-egger
    output3<-paste0(outputpath, outcomefilename,'_',metaname,'_QC_outlierSNPs') #output data outlier snp
    
    dat <- harmonise_data(exposure_dat = exposure_clumped,outcome_dat = outcome_dat_0,action=2)
    data3<-format_radial(dat$beta.exposure,dat$beta.outcome,dat$se.exposure,dat$se.outcome,dat$SNP)
    res1 <- tryCatch(
            expr = {
            ivw_radial(data3, 0.05, 1, 0.0001)
            },
            error = function(e) {
            # 如果出现错误，打印错误信息并跳过当前循环
            message("ivw_radial make mistaks,continue...")
            return(NULL)
            }
        )

    if (!is.null(res1)) {
    res2<-egger_radial(data3,0.05,1)
    write.table(as.data.frame(res1$data),output1,quote=F,row.names=F,col.names=T,sep="\t")
    write.table(as.data.frame(res2$data),output2,quote=F,row.names=F,col.names=T,sep="\t")
    if (class(res1$outliers)=='data.frame'){
        a1<-res1$outliers$SNP
    }else{
        a1<-c('None')
    }
    if (class(res2$outliers)=='data.frame'){
        a2<-res2$outliers$SNP
    }else{
        a2<-c('None')
    }
    a3<-unique(c(a1,a2))
    write.table(data.frame(a3),output3,quote=F,row.names=F,col.names=F,sep="\t")
    
      ######################################## 4. sensitivity analysis #########################################
      output1<-paste0(outputpath, outcome,'_',metaname,'_dat')
      write.table(dat, output1, quote=F, row.names=F,col.names=T,sep="\t")
      library(MRPRESSO)
      output1<-paste0(outputpath, outcome,'_',metaname,'_sensitivity_MR_PRESSO_distortion') #distortion_test
      output2<-paste0(outputpath, outcome,'_',metaname,'_sensitivity_MR_PRESSO_global') #global_test
      output3<-paste0(outputpath, outcome,'_',metaname,'_sensitivity_MR_PRESSO_main') #main_MR_results
      output4<-paste0(outputpath, outcome,'_',metaname,'_sensitivity_MR_PRESSO_outlierTest') #eachsnp.outlier_test
      bzx<-dat$beta.exposure
      sebzx<-dat$se.exposure
      bzy<-dat$beta.outcome
      sebzy<-dat$se.outcome
      data3<-as.data.frame(cbind(bzx,sebzx,bzy,sebzy))
      results<-mr_presso(BetaOutcome = "bzy", BetaExposure = "bzx", SdOutcome= "sebzy", SdExposure = "sebzx", OUTLIERtest = TRUE, DISTORTIONtest= TRUE, data = data3, NbDistribution = 10000, SignifThreshold = 0.05)
      a1<-results$`MR-PRESSO results`$`Distortion Test` #testing of significant distortion in the causal estimates before and after outlier removal
      a2<-results$`Main MR results`
      a3<-results$`MR-PRESSO results`$`Global Test`$Pvalue #detection of horizontal pleiotropy
      a4<-results$`MR-PRESSO results`$`Global Test`$RSSobs
      a5<-results$`MR-PRESSO results`$`Outlier Test` #correction of horizontal pleiotropy via outlier removal
      b2<-cbind(a4,a3)
      a2A<-cbind(a2[2],a2[3],a2[4],a2[5],a2[6])
      write.table(a1,output1,quote=F,row.names=F,col.names=F,sep="\t")
      write.table(b2,output2,quote=F,row.names=F,col.names=F,sep="\t")
      write.table(a2A,output3,quote=F,row.names=T,col.names=T,sep="\t")
      write.table(a5,output4,quote=F,row.names=F,col.names=T,sep="\t")
  
      output1 <- paste0(outputpath, outcome,'_',metaname,'_sensitivity_MR_egger_intercept')
      Inter <- mr_egger_regression(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
      out <- cbind(Inter$b_i, Inter$se_i, Inter$pval_i) #egger intercept
      write.table(out,output1,quote=F,row.names=F,col.names=F,sep="\t") #egger intercept
  
  
      output1 <- paste0(outputpath, outcome,'_',metaname,'_sensitivity_MR_LOO')
      loo <- mr_leaveoneout(dat)
      b <- data.frame(loo$SNP,loo$b,loo$se,loo$p)
      write.table(b,output1,quote=F,row.names=F,col.names=T,sep="\t") #leave one out results
  
      output1 <- paste0(outputpath, outcome,'_',metaname,'_sensitivity_MR_steigerDirection')
      out <- directionality_test(dat)
      write.table(out, output1,quote=F,row.names=F,col.names=T,sep="\t") #steigerDirection
      
  }
}
}