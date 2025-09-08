library(TwoSampleMR)
library(data.table)
library('ieugwasr')
library(RadialMR)

use_data <- read.csv("BONsig_bing.csv")

meta_list <- unique(use_data$meta)
diseases_list <- unique(use_data$diseases)

meta_list_use <- meta_list
diseases_list_use <- diseases_list

bfile='EUR'

for(i in seq(1,length(meta_list_use))){
  meta = meta_list_use[i]
  print("**************************")
  print(i)
  print("**************************")
  
  threshold=5e-8
  threshold_folder="MR_5e8"
  
  exposure_dat_0 <- read_exposure_data(
    filename = paste0(meta,'/GWAS_merge.glm.linear'),
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
  exposure_dat_0$exposure <- rep(meta,nrow(exposure_dat_0))
  # https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  exposure_dat_0$rsid <- exposure_dat_0$SNP
  exposure_dat_0$trait_id <- rep(meta,nrow(exposure_dat_0))
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
  list_choose = use_data[use_data$meta==meta,]
  diseases_choose = list_choose$diseases
  
  
for(j in seq(1,length(diseases_choose))){
    print(j)
    out = diseases_choose[j]

    outputpath=paste0(meta,'/',out,'/')
    filepath = paste0(outputpath, out, "_",meta,"_MR_IVW")

    if (dir.exists(outputpath)) {
        # 如果存在，则跳过当前循环
        message(paste("Directory already exists:", outputpath))
        if(file.exists(filepath)){
        next
        }else{
        unlink(outputpath, recursive = TRUE)
        }
        
    } else {
        outcomefile <- paste0(out)
        
        if (file.exists(outcomefile)) {
        # 如果不存在，则创建目录
        dir.create(outputpath, recursive=TRUE)
          snps=exposure_clumped$SNP
          outcome_dat_0 <- read_outcome_data(
            snps = snps,
            filename = outcomefile,
            sep = "\t",
            snp_col = "rsids",
            beta_col = "beta",
            se_col = "sebeta",
            effect_allele_col = "alt",
            other_allele_col = "ref",
            eaf_col = "af_alt",  # merge列
            pval_col = "mlogp",
            log_pval = TRUE
          )

        
        outcome = as.character(out)
        outcome_dat_0$outcome <- rep(outcome,nrow(outcome_dat_0))
        sampleSize_outcome=243096
        outcome_dat_0$samplesize.outcome <- rep(sampleSize_outcome,nrow(outcome_dat_0))

        outcomefilename = outcome
        output1<-paste0(outputpath, outcomefilename,'_',meta,'_QC_IVW') #output data IVW
        output2<-paste0(outputpath, outcomefilename,'_',meta,'_QC_MR-egger') #output data MR-egger
        output3<-paste0(outputpath, outcomefilename,'_',meta,'_QC_outlierSNPs') #output data outlier snp
        
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
            
            
            ######################################### 3. MR analysis #########################################
            output1 <- paste0(outputpath, outcome,'_',meta,'_MR_WM') #weight median
            output2 <- paste0(outputpath, outcome,'_',meta,'_MR_IVW') #ivw
            output3 <- paste0(outputpath, outcome,'_',meta,'_MR_egger') #egger slope
            output4 <- paste0(outputpath, outcome,'_',meta,'_MR_W') #weight mode
            output5 <- paste0(outputpath, outcome,'_',meta,'_MR_robustAdjProfile') #Robust adjusted profile score
            output6 <- paste0(outputpath, outcome,'_',meta,'_MR_waldRatio') #Wald Ratio
            
            
            # remove outlier snps after quality control
            dat <- dat[!(dat$SNP %in% a3),]
            source("/public/share/tmp/mr_modified.R")
            tsmr1<-mr(dat, method_list=c("mr_weighted_median"))
            tsmr2<-mr(dat, method_list=c("mr_ivw"))
            tsmr5<-mr_modified(dat, method_list=c("mr_raps"))
            tsmr6<-mr(dat, method_list=c("mr_wald_ratio"))
            tsmr3<-mr(dat, method_list=c("mr_egger_regression"))
            tsmr4<-mr(dat, method_list=c("mr_weighted_mode"))
            a1<-cbind(outcome, meta,tsmr1$nsnp,tsmr1$b,(tsmr1$b-1.96*tsmr1$se),(tsmr1$b+1.96*tsmr1$se),tsmr1$pval) #weight median
            a2<-cbind(outcome, meta,tsmr2$nsnp,tsmr2$b,(tsmr2$b-1.96*tsmr2$se),(tsmr2$b+1.96*tsmr2$se),tsmr2$pval) #ivw
            a5<-cbind(outcome, meta,tsmr5$nsnp,tsmr5$b,(tsmr5$b-1.96*tsmr5$se),(tsmr5$b+1.96*tsmr5$se),tsmr5$pval) #raps
            a6<-cbind(outcome, meta,tsmr6$nsnp,tsmr6$b,(tsmr6$b-1.96*tsmr6$se),(tsmr6$b+1.96*tsmr6$se),tsmr6$pval) #WR
            CI <- 0.95
            lowerCI <- function(beta,df,SE){
            return(beta - (qt((1-CI)/2, df, lower.tail = FALSE) * SE))
            }
            upperCI <- function(beta,df,SE){
            return(beta + (qt((1-CI)/2, df, lower.tail = FALSE) * SE))
            }
            a3 <- cbind(outcome, meta,tsmr3$nsnp,tsmr3$b, mapply(lowerCI, tsmr3$b, tsmr3$nsnp - 2, tsmr3$se), mapply(upperCI, tsmr3$b, tsmr3$nsnp - 2, tsmr3$se), tsmr3$pval) #egger slope
            a4 <- cbind(outcome, meta,tsmr4$nsnp,tsmr4$b, mapply(lowerCI, tsmr4$b, tsmr4$nsnp - 1, tsmr4$se), mapply(upperCI, tsmr4$b, tsmr4$nsnp - 1, tsmr4$se), tsmr4$pval) #weight mode
            write.table(a1,output1,quote=F,row.names=F,col.names=F,sep="\t") #weight median
            write.table(a2,output2,quote=F,row.names=F,col.names=F,sep="\t") #ivw
            write.table(a5,output5,quote=F,row.names=F,col.names=F,sep="\t") #raps
            write.table(a6,output6,quote=F,row.names=F,col.names=F,sep="\t") #wr
            write.table(a3,output3,quote=F,row.names=F,col.names=F,sep="\t") #egger slope
            write.table(a4,output4,quote=F,row.names=F,col.names=F,sep="\t") #weight mode
        } else{
            unlink(outputpath, recursive = TRUE)
            next
        }
        
        }else {
        message(paste("Disease does not exist:", out))
        next
        }
    }
}
}