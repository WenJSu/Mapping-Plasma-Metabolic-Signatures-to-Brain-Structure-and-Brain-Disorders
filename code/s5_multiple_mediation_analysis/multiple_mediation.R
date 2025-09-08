library(lavaan)
library(semTools)
library(dplyr)
library(data.table)

fill_with_median <- function(x) {
  if (any(is.na(x))) {
    median_value <- median(x, na.rm = TRUE)
    x[is.na(x)] <- median_value
  }
  return(x)
}

imagedata = read.csv('imgdata_aft_pro.csv')
eid = imagedata['eid']

metadata <- fread("meta_aft_pro.csv", sep = ",", header = TRUE)
metadata_merge <- merge(eid,metadata,by='eid')

covadata <- fread("cova_data_more_withPC_all2visits.csv", sep = ",", header = TRUE)
covadata_used <- covadata[,c('eid','age2','Sex')]

cov_all <- covadata_used
columns_to_standardize <- c('age2')
temp <- as.data.frame(scale(cov_all[,..columns_to_standardize]))
colnames(temp)<- c('age2')
cov_all[,columns_to_standardize] = temp[,columns_to_standardize]

refer_list <- fread("combined_data2_readj_sig.csv", sep = ",", header = TRUE)
usedata_or <- merge(imagedata, metadata_merge, by=c('eid'))
usedata_or <- merge(usedata_or, cov_all, by=c('eid'))

targetfiles <- unique(refer_list$diseases)


target_path <- "/data/"
csv_paths <- list.files(
  path = target_path,
  pattern = "\\.csv$",
  full.names = FALSE  # 不返回完整路径，仅返回csv文件名称
)
csv_list <- csv_paths

for(FID in seq(1,length(csv_list))){
  print(paste0("FID:",FID))
  FID_use = csv_list[FID]
  print(FID_use)
  FID_dir <- paste0(target_path,FID_use)
  FID_df <- read.csv(FID_dir)
  colnames(FID_df)<-c("eid","lifestyle")
  FID_df$lifestyle <- scale(FID_df$lifestyle)
  for(i in seq(1,length(targetfiles))){
  
    num_list_ill <- c()
    num_XNA_list <- c()
    CFI <- c()
    RMSEA <- c()
    SRMR <- c()
  
    list_disease <- c()
    list_meta <- c()
    list_IDP <- c()
  
    list_a <- c()
    list_b <- c()
    list_ab <- c()
    list_cp <- c()
    list_d <- c()
    list_ep <- c()
    list_abd <- c()
    list_total <- c()
    list_dtotal <- c()
    list_all <- c()
    
    list_ap <- c()
    list_bp <- c()
    list_abp <- c()
    list_cpp <- c()
    list_dp <- c()
    list_epp <- c()
    list_abdp <- c()
    list_totalp <- c()
    list_dtotalp <- c()
    list_allp <- c()
    
    fID <- c()
  
    Padj_FDR_a <- c()
    Padj_FDR_b <- c()
    Padj_FDR_c <- c()
    Padj_FDR_ab <-c()
    Padj_FDR_total <- c()
  
    list_IDPcate <- c()
    list_IDPcate_sum <- c()
  
    print(i)
    ill_use = targetfiles[i]
    print(ill_use)
    use_ill_meta <- refer_list[refer_list$diseases == ill_use, ]
  
    targetdata <- read.csv(paste0(ill_use,'.csv'),sep=',',header=T)
    targetdata <- targetdata[,c('eid','target_y')]
    targetdata = na.omit(targetdata)
  
    usedata <- merge(usedata_or, targetdata, by=c('eid'))
  
    
    for(t in seq(1,nrow(use_ill_meta))){
      run_process <- paste0(t, "/", nrow(use_ill_meta))
      print(run_process)
    
  
      list_disease <- c(list_disease, ill_use)
      list_meta <- c(list_meta, use_ill_meta$meta[t])
      list_IDP <- c(list_IDP, use_ill_meta$IDP[t])
      fID <- c(fID,FID_use)
  
      IDP_unique = use_ill_meta$IDP[t]
      IDP_unique <- strsplit(IDP_unique, ",")[[1]]
      
      cols_to_extract <- c("eid", use_ill_meta$meta[t], IDP_unique, "age2", "Sex","target_y" )
      usedata_choose <- usedata[, cols_to_extract]
      colnames(usedata_choose)[colnames(usedata_choose) == use_ill_meta$meta[t]] <- "meta"
      
      usedata_choose <- merge(usedata_choose, FID_df, by=c('eid'))
      

      usedata_choose_XNA <- usedata_choose %>%
        mutate_all(fill_with_median)
      # usedata_choose_XNA <- na.omit(usedata_choose)
      num_XNA_list <- c(num_XNA_list, nrow(usedata_choose_XNA))
      num_list_ill <- c(num_list_ill, sum(usedata_choose_XNA$target_y == 1, na.rm = TRUE))
      if(length(IDP_unique)==1){
        usedata_choose_XNA$IDP_mean <- usedata_choose_XNA[[IDP_unique]]
      }else{
        usedata_choose_XNA$IDP_mean <- rowMeans(usedata_choose_XNA[, IDP_unique], na.rm = TRUE)
      }
      
      model <- paste(
        "meta ~ d*lifestyle+age2+Sex",
        "IDP_mean ~ a*meta + age2 + Sex",  # meta 影响 IDP (路径 a)
        "target_y ~ b*IDP_mean + age2 + Sex",  # meta 和 IDP 同时预测 ill
        "target_y ~ cp*meta + age2 + Sex",
        "target_y ~ ep*lifestyle + age2 + Sex",
        "ab := a * b",
        "abd := a*b*d",
        "total := cp + ab",
        "dtotal := d*total",
        "all:=ep+dtotal",
        sep = '\n'
      )
      sem_result <- sem(
        model,
        data = usedata_choose_XNA,
        ordered = c("target_y"),    # 确保变量名用引号包裹，且存在于数据中
        estimator = "DWLS",       # 或 estimator = "DWLS"
        se = "bootstrap",
        bootstrap = 1000
      )
      result <- summary(sem_result, fit.measures = TRUE, standardized=TRUE)
      
      CFI <- c(CFI, result[["fit"]][["cfi"]])
      RMSEA <- c(RMSEA, result[["fit"]][["rmsea"]])
      SRMR <- c(SRMR, result[["fit"]][["srmr"]])
      
      list_a <- c(list_a, result[["pe"]][["est"]][4])
      list_b <- c(list_b, result[["pe"]][["est"]][7])
      list_ab <- c(list_ab, result[["pe"]][["est"]][27])
      list_cp <- c(list_cp, result[["pe"]][["est"]][10])
      list_d <- c(list_d, result[["pe"]][["est"]][1])
      list_ep <- c(list_ep, result[["pe"]][["est"]][11])
      list_abd <- c(list_abd, result[["pe"]][["est"]][28])
      list_total <- c(list_total, result[["pe"]][["est"]][29])
      list_dtotal <- c(list_dtotal, result[["pe"]][["est"]][30])
      list_all <- c(list_all, result[["pe"]][["est"]][31])
      
      list_ap <- c(list_ap, result[["pe"]][["pvalue"]][4])
      list_bp <- c(list_bp, result[["pe"]][["pvalue"]][7])
      list_abp <- c(list_abp, result[["pe"]][["pvalue"]][27])
      list_cpp <- c(list_cpp, result[["pe"]][["pvalue"]][10])
      list_dp <- c(list_dp, result[["pe"]][["pvalue"]][1])
      list_epp <- c(list_epp, result[["pe"]][["pvalue"]][11])
      list_abdp <- c(list_abdp, result[["pe"]][["pvalue"]][28])
      list_totalp <- c(list_totalp, result[["pe"]][["pvalue"]][29])
      list_dtotalp <- c(list_dtotalp, result[["pe"]][["pvalue"]][30])
      list_allp <- c(list_allp, result[["pe"]][["pvalue"]][31])
    }
  
  
  resultframe <- data.frame(diseases=list_disease, meta=list_meta, IDP=list_IDP,
                            num_XNA=num_XNA_list, num_XNA_ill=num_list_ill, CFI=CFI, RMSEA=RMSEA, SRMR=SRMR,
                            a=list_a, a_pvalues=list_ap, b=list_b, b_pvalues=list_bp, 
                            ab=list_ab, ab_pvalues=list_abp,cp=list_cp, cp_pvalues=list_cpp, 
                            d=list_d, d_pvalues=list_dp,ep=list_ep, ep_pvalues=list_epp, 
                            abd=list_abd, abd_pvalues=list_abdp,total=list_total, total_pvalues=list_totalp, 
                            dtotal=list_dtotal, dtotal_pvalues=list_dtotalp, 
                            all=list_all, all_pvalues=list_allp, fID=fID)
  
  new_dir <- paste0(FID_use,"/")
  dir.create(new_dir, recursive=TRUE)
  output_dir <- paste0(new_dir, ill_use,".csv")
  write.table(resultframe, output_dir,sep=',',row.names = F)
  }
}

print("Finish!!!")
