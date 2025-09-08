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
cov_all[,columns_to_standardize] = temp[,columns_to_standardize]

refer_list <- fread("pair_bing.csv", sep = ",", header = TRUE)
usedata_or <- merge(imagedata, metadata_merge, by=c('eid'))
usedata_or <- merge(usedata_or, cov_all, by=c('eid'))

targetfiles <- unique(refer_list$diseases)
targetfiles <- targetfiles

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
  list_c <- c()
  list_ab <- c()
  list_total <- c()

  list_apvalues <- c()
  list_bpvalues <- c()
  list_cpvalues <- c()
  list_abpvalues <- c()
  list_totalpvalues <- c()

  Padj_FDR_a <- c()
  Padj_FDR_b <- c()
  Padj_FDR_c <- c()
  Padj_FDR_ab <-c()
  Padj_FDR_total <- c()

  list_IDPcate <- c()
  list_IDPcate_sum <- c()

  print(i)
  ill_use = targetfiles[i]
  use_ill_meta <- refer_list[refer_list$diseases == ill_use, ]

  targetdata <- read.csv(paste0(ill_use,'.csv'),sep=',',header=T)
  targetdata <- targetdata[,c('eid','target_y')]
  targetdata = na.omit(targetdata)

  usedata <- merge(usedata_or, targetdata, by=c('eid'))
  meta_list = unique(use_ill_meta$meta)

  
  for(t in seq(1,length(meta_list))){
    run_process <- paste0(t, "/", length(meta_list))
    print(run_process)
    
    meta_use = meta_list[t]
    use_ill_meta_meta <- use_ill_meta[use_ill_meta$meta == meta_use, ]
    


    print("cortical thickness")
    use_cate_df_dire_IDP <- use_ill_meta_meta %>%
      filter(Category.name %in% c("cortical thickness"))
    if(nrow(use_cate_df_dire_IDP) > 0){
      print("in cortical thickness")

      use_cate_df_dire_IDP$IDP = paste0("X", use_cate_df_dire_IDP$UKBID.x, ".2.0")
      list_disease <- c(list_disease, ill_use)
      list_meta <- c(list_meta, meta_use)
      
      IDP_unique <- unique(use_cate_df_dire_IDP$IDP)
      IDPcate_unique <- unique(use_cate_df_dire_IDP$Category.name)
      

      list_IDP <- c(list_IDP, paste(IDP_unique, collapse = ","))
      list_IDPcate <- c(list_IDPcate, paste(IDPcate_unique, collapse = ","))
      list_IDPcate_sum <- c(list_IDPcate_sum, "cortical thickness")
      
      cols_to_extract <- c("eid", meta_use, IDP_unique, "age2", "Sex","target_y" )
      usedata_choose <- usedata[, cols_to_extract]
      colnames(usedata_choose)[colnames(usedata_choose) == meta_use] <- "meta"
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
        "IDP_mean ~ a*meta + age2 + Sex",  # meta 影响 IDP (路径 a)
        "target_y ~ b*IDP_mean + age2 + Sex",  # meta 和 IDP 同时预测 ill
        "target_y ~ cp*meta + age2 + Sex",
        "ab := a * b",
        "total := cp + ab", 
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
      
      list_a <- c(list_a, result[["pe"]][["est"]][1])
      list_b <- c(list_b, result[["pe"]][["est"]][4])
      list_c <- c(list_c, result[["pe"]][["est"]][7])
      list_ab <- c(list_ab, result[["pe"]][["est"]][21])
      list_total <- c(list_total, result[["pe"]][["est"]][22])
      
      list_apvalues <- c(list_apvalues, result[["pe"]][["pvalue"]][1])
      list_bpvalues <- c(list_bpvalues, result[["pe"]][["pvalue"]][4])
      list_cpvalues <- c(list_cpvalues, result[["pe"]][["pvalue"]][7])
      list_abpvalues <- c(list_abpvalues, result[["pe"]][["pvalue"]][21])
      list_totalpvalues <- c(list_totalpvalues, result[["pe"]][["pvalue"]][22])
    }
    
    #2. cortical area
    print("cortical area")
    use_cate_df_dire_IDP <- use_ill_meta_meta %>%
      filter(Category.name %in% c("cortical area"))
    if(nrow(use_cate_df_dire_IDP) > 0){
      print("in cortical area")
      use_cate_df_dire_IDP$IDP = paste0("X", use_cate_df_dire_IDP$UKBID.x, ".2.0")
      list_disease <- c(list_disease, ill_use)
      list_meta <- c(list_meta, meta_use)
      
      IDP_unique <- unique(use_cate_df_dire_IDP$IDP)
      IDPcate_unique <- unique(use_cate_df_dire_IDP$Category.name)
      
      
      list_IDP <- c(list_IDP, paste(IDP_unique, collapse = ","))
      list_IDPcate <- c(list_IDPcate, paste(IDPcate_unique, collapse = ","))
      list_IDPcate_sum <- c(list_IDPcate_sum, "cortical thickness")
      
      cols_to_extract <- c("eid", meta_use, IDP_unique, "age2", "Sex","target_y" )
      usedata_choose <- usedata[, cols_to_extract]
      colnames(usedata_choose)[colnames(usedata_choose) == meta_use] <- "meta"
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
        "IDP_mean ~ a*meta + age2 + Sex",  # meta 影响 IDP (路径 a)
        "target_y ~ b*IDP_mean + age2 + Sex",  # meta 和 IDP 同时预测 ill
        "target_y ~ cp*meta + age2 + Sex",
        "ab := a * b",
        "total := cp + ab", 
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
      
      list_a <- c(list_a, result[["pe"]][["est"]][1])
      list_b <- c(list_b, result[["pe"]][["est"]][4])
      list_c <- c(list_c, result[["pe"]][["est"]][7])
      list_ab <- c(list_ab, result[["pe"]][["est"]][21])
      list_total <- c(list_total, result[["pe"]][["est"]][22])
      
      list_apvalues <- c(list_apvalues, result[["pe"]][["pvalue"]][1])
      list_bpvalues <- c(list_bpvalues, result[["pe"]][["pvalue"]][4])
      list_cpvalues <- c(list_cpvalues, result[["pe"]][["pvalue"]][7])
      list_abpvalues <- c(list_abpvalues, result[["pe"]][["pvalue"]][21])
      list_totalpvalues <- c(list_totalpvalues, result[["pe"]][["pvalue"]][22])
    }
    #3.regional and tissue volume
    print("regional and tissue volume")
    use_cate_df_dire_IDP <- use_ill_meta_meta %>%
      filter(Category.name %in% c("regional and tissue volume"))
    if(nrow(use_cate_df_dire_IDP) > 0){
      print("in regional and tissue volume")
      use_cate_df_dire_IDP$IDP = paste0("X", use_cate_df_dire_IDP$UKBID.x, ".2.0")
      list_disease <- c(list_disease, ill_use)
      list_meta <- c(list_meta, meta_use)
      
      IDP_unique <- unique(use_cate_df_dire_IDP$IDP)
      IDPcate_unique <- unique(use_cate_df_dire_IDP$Category.name)
      
      
      list_IDP <- c(list_IDP, paste(IDP_unique, collapse = ","))
      list_IDPcate <- c(list_IDPcate, paste(IDPcate_unique, collapse = ","))
      list_IDPcate_sum <- c(list_IDPcate_sum, "cortical thickness")
      
      cols_to_extract <- c("eid", meta_use, IDP_unique, "age2", "Sex","target_y" )
      usedata_choose <- usedata[, cols_to_extract]
      colnames(usedata_choose)[colnames(usedata_choose) == meta_use] <- "meta"
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
        "IDP_mean ~ a*meta + age2 + Sex",  # meta 影响 IDP (路径 a)
        "target_y ~ b*IDP_mean + age2 + Sex",  # meta 和 IDP 同时预测 ill
        "target_y ~ cp*meta + age2 + Sex",
        "ab := a * b",
        "total := cp + ab", 
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
      
      list_a <- c(list_a, result[["pe"]][["est"]][1])
      list_b <- c(list_b, result[["pe"]][["est"]][4])
      list_c <- c(list_c, result[["pe"]][["est"]][7])
      list_ab <- c(list_ab, result[["pe"]][["est"]][21])
      list_total <- c(list_total, result[["pe"]][["est"]][22])
      
      list_apvalues <- c(list_apvalues, result[["pe"]][["pvalue"]][1])
      list_bpvalues <- c(list_bpvalues, result[["pe"]][["pvalue"]][4])
      list_cpvalues <- c(list_cpvalues, result[["pe"]][["pvalue"]][7])
      list_abpvalues <- c(list_abpvalues, result[["pe"]][["pvalue"]][21])
      list_totalpvalues <- c(list_totalpvalues, result[["pe"]][["pvalue"]][22])
    }
    #3.FA
    print("WM tract FA")
    use_cate_df_dire_IDP <- use_ill_meta_meta %>%
      filter(Category.name %in% c("WM tract FA"))
    if(nrow(use_cate_df_dire_IDP) > 0){
      print("in WM tract FA")
      use_cate_df_dire_IDP$IDP = paste0("X", use_cate_df_dire_IDP$UKBID.x, ".2.0")
      list_disease <- c(list_disease, ill_use)
      list_meta <- c(list_meta, meta_use)
      
      IDP_unique <- unique(use_cate_df_dire_IDP$IDP)
      IDPcate_unique <- unique(use_cate_df_dire_IDP$Category.name)
      
      
      list_IDP <- c(list_IDP, paste(IDP_unique, collapse = ","))
      list_IDPcate <- c(list_IDPcate, paste(IDPcate_unique, collapse = ","))
      list_IDPcate_sum <- c(list_IDPcate_sum, "cortical thickness")
      
      cols_to_extract <- c("eid", meta_use, IDP_unique, "age2", "Sex","target_y" )
      usedata_choose <- usedata[, cols_to_extract]
      colnames(usedata_choose)[colnames(usedata_choose) == meta_use] <- "meta"
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
        "IDP_mean ~ a*meta + age2 + Sex",  # meta 影响 IDP (路径 a)
        "target_y ~ b*IDP_mean + age2 + Sex",  # meta 和 IDP 同时预测 ill
        "target_y ~ cp*meta + age2 + Sex",
        "ab := a * b",
        "total := cp + ab", 
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
      
      list_a <- c(list_a, result[["pe"]][["est"]][1])
      list_b <- c(list_b, result[["pe"]][["est"]][4])
      list_c <- c(list_c, result[["pe"]][["est"]][7])
      list_ab <- c(list_ab, result[["pe"]][["est"]][21])
      list_total <- c(list_total, result[["pe"]][["est"]][22])
      
      list_apvalues <- c(list_apvalues, result[["pe"]][["pvalue"]][1])
      list_bpvalues <- c(list_bpvalues, result[["pe"]][["pvalue"]][4])
      list_cpvalues <- c(list_cpvalues, result[["pe"]][["pvalue"]][7])
      list_abpvalues <- c(list_abpvalues, result[["pe"]][["pvalue"]][21])
      list_totalpvalues <- c(list_totalpvalues, result[["pe"]][["pvalue"]][22])
    }
  }
Padj_FDR_a <- c(Padj_FDR_a, p.adjust(list_apvalues,method='BH'))
Padj_FDR_b <- c(Padj_FDR_b, p.adjust(list_bpvalues,method='BH'))
Padj_FDR_c <- c(Padj_FDR_c, p.adjust(list_cpvalues,method='BH'))
Padj_FDR_ab <-c(Padj_FDR_ab, p.adjust(list_abpvalues,method='BH'))
Padj_FDR_total <- c(Padj_FDR_total, p.adjust(list_totalpvalues,method='BH'))

resultframe <- data.frame(diseases=list_disease, meta=list_meta, IDP=list_IDP,
                          num_XNA=num_XNA_list, num_XNA_ill=num_list_ill, CFI=CFI, RMSEA=RMSEA, SRMR=SRMR, a_estimate=list_a,
                          a_pvalues=list_apvalues, a_Padj_FDR=Padj_FDR_a, b_estimate=list_b,
                          b_pvalues=list_bpvalues, b_Padj_FDR=Padj_FDR_b, c_estimate=list_c,
                          c_pvalues=list_cpvalues, c_Padj_FDR=Padj_FDR_c,
                          ab_estimate=list_ab,
                          ab_pvalues=list_abpvalues, ab_Padj_FDR=Padj_FDR_ab, total_estimate=list_total,
                          total_pvalues=list_totalpvalues, total_Padj_FDR=Padj_FDR_total,
                          mean_IDP_cate=list_IDPcate, mean_IDP_cate_sum=list_IDPcate_sum)

output_dir <- paste0(ill_use,".csv")
write.table(resultframe, output_dir,sep=',',row.names = F)
}

print("Finish!!!")
