# survial analysis
################################## with more covariates
library(survival)
library("fastDummies")
library(data.table)
library(dplyr)
library(tidyr)

metadata <- fread("meta_aft_pro.csv", sep = ",", header = TRUE)
meta_index <- colnames(metadata)


diseasedata <- fread("diseases.csv", sep = ",", header = TRUE)
disease_index <- diseasedata$diseases

cov2 <- fread("Covariates.csv", sep = ",", header = TRUE)
cov2 <- cov2[,c('eid','Race','FastingTime','Statin')]
cov2_used <- fastDummies::dummy_cols(cov2, select_columns = "Race", remove_first_dummy = TRUE)

covadata <- fread("cova_data_more_withPC_all2visits.csv", sep = ",", header = TRUE)
covadata_used <- covadata[,c('eid','age0','Sex','Townsend_index','smoking0','drinking0','BMI0','Education')]
covadata_used <- covadata_used[covadata_used$smoking0!=-3 & covadata_used$drinking0!=-3,]
covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "drinking0", remove_first_dummy = TRUE)
covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "smoking0", remove_first_dummy = TRUE)

cov_all <- merge(cov2_used,covadata_used,by='eid')
columns_to_standardize <- c('age0', 'Townsend_index', 'BMI0', 'Education','FastingTime')

temp <- as.data.frame(scale(cov_all[,..columns_to_standardize]))
cov_all[,columns_to_standardize] = temp[,columns_to_standardize]

# cox harzard model, with delta as continuous variable, in all subjs, with more covariates
list_meta <- c()
list_targetfile <- c()
list_numberAll <- c()
list_numberTarget <- c()
list_pvalues <- c()
list_zvalues <- c()
list_coef <- c()
list_expcoef <- c()
list_secoef <- c()
list_FDR_separate <- c()
list_up <- c()
list_down <- c()
for(i in seq(2,length(meta_index))){
  out1 <- paste0(i-1,"/",length(meta_index)-1)
  print("******************************")
  print(out1)
  print("******************************")
  meta_name <- meta_index[i]
  meta_single <- subset(metadata, select = c(1, i))
  pvalues <- c()


  for(j in seq(1,length(disease_index))){
    if(j%%10==0){
      out2 <- paste0(j,"/",length(disease_index))
      print(out2)
    }
    targetfile <- disease_index[j]
    list_meta <- c(list_meta, meta_index[i])
    list_targetfile <- c(list_targetfile, targetfile)
    disease_path <- paste0(targetfile,".csv")
    disease_single <- fread(disease_path, sep = ",", header = TRUE)
    disease_single <- disease_single[,c("eid", "target_y", "BL2Target_yrs")]
    
    # convert variables merge
    tempdata <- merge(meta_single, disease_single, by=c('eid'))
    useddata <- merge(tempdata,cov_all, by=c('eid'))
    
    useddata$time <- useddata$BL2Target_yrs
    useddata <- useddata[useddata$time>0,]
    useddata=na.omit(useddata)
    
    list_numberAll <- c(list_numberAll, nrow(useddata))
    list_numberTarget <- c(list_numberTarget, nrow(useddata[useddata$target_y==1,]))
    
    # cross-section analysis
    cleaned_data <- useddata %>% select(meta_name,"BL2Target_yrs","target_y","Race_2","Race_3","Race_4","age0","Sex","Townsend_index","BMI0","Education","drinking0_1","drinking0_2","smoking0_1","smoking0_2",'FastingTime','Statin')
    colnames(cleaned_data) <- c("meta","years","disease","Race_2","Race_3","Race_4","age0","Sex","TDI","BMI0","Education","drinking0_1","drinking0_2","smoking0_1","smoking0_2",'FastingTime','Statin')
    
    cox_fit <- coxph(Surv(years, disease) ~ meta+Race_2+Race_3+Race_4+age0+Sex+TDI+BMI0+Education+smoking0_1+smoking0_2+drinking0_1+drinking0_2+FastingTime+Statin, data=cleaned_data)

    result <- summary(cox_fit)$coefficients
    CI_range <- summary(cox_fit)$conf.int[1,c(3,4)]
    # Check for the possibility of infinite coefficients
    if (any(is.infinite(coefficients(cox_fit)))) {
      print(j)
      print("Warning: Model coefficients may be infinite.")
    }
    list_pvalues <- c(list_pvalues, result[1,5])
    list_zvalues <- c(list_zvalues, result[1,4])
    list_coef <- c(list_coef, result[1,1])
    list_expcoef <- c(list_expcoef, result[1,2])
    list_secoef <- c(list_secoef, result[1,3])
    list_up <- c(list_up, CI_range[2])
    list_down <- c(list_down, CI_range[1])
    pvalues <- c(pvalues, result[1,5])
  }
  print(i)
  list_FDR_separate <- c(list_FDR_separate,p.adjust(pvalues, method = 'BH'))
}

resultframe <- data.frame(meta=list_meta, target=list_targetfile, numberAll=list_numberAll, 
                          numberTarget=list_numberTarget, expcoef=list_expcoef, coef=list_coef, 
                          down=unname(list_down), up=unname(list_up), secoef=list_secoef, zvalue=list_zvalues, 
                          pvalue=list_pvalues, separateFDR=list_FDR_separate)

resultframe$Padj_FDR_overall <- p.adjust(resultframe$pvalue,method='BH')
write.table(resultframe, 'cox_meta_diseases.csv',sep=',',row.names = F)
print("Finish!")