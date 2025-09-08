library(data.table)
library(dplyr)
library(tidyr)

# load meta and img data
metadata <- fread("meta_aft_pro.csv", sep = ",", header = TRUE)
imagedata <- fread('img_pro_ins2_cov.csv', sep = ",", header = TRUE)

meta_index <- colnames(metadata)
image_index <- colnames(imagedata)[1:269]

# load excluding_df
ex_df <- fread("excludingSubjs.csv", sep = ",", header = TRUE)

# load covariates
cov2 <- fread("Covariates.csv", sep = ",", header = TRUE)
cov2 <- cov2[,c('eid','Race','FastingTime','Statin')]
cov2_used <- fastDummies::dummy_cols(cov2, select_columns = "Race", remove_first_dummy = TRUE)

allSiteInfo <- read.table('eTIV_site_new2.csv', sep=',' ,header=T, check.names=F)
allSiteInfo <- allSiteInfo[,c('eid','p26521_i2','Site')]
allSiteInfo <- allSiteInfo[allSiteInfo$'Site'!='11028',] # remove the site with too little subjects
allSiteInfo <- na.omit(allSiteInfo)
colnames(allSiteInfo)[which(colnames(allSiteInfo) == "p26521_i2")] <- "eTIV"
allSiteInfo$Site <- as.character(allSiteInfo$Site)
allSiteInfo <- fastDummies::dummy_cols(allSiteInfo, select_columns = "Site", remove_first_dummy = TRUE)


imagedata <- imagedata[imagedata$smoking2!=-3 & imagedata$drinking2!=-3,]
imagedata <- fastDummies::dummy_cols(imagedata, select_columns = "drinking2", remove_first_dummy = TRUE)
imagedata <- fastDummies::dummy_cols(imagedata, select_columns = "smoking2", remove_first_dummy = TRUE)

cov_all <- merge(cov2_used,allSiteInfo,by='eid')
cov_inner = imagedata[,c("eid","Sex","Townsend_index","Education","BMI2","drinking2_1","drinking2_2","smoking2_1","smoking2_2","p21003_i2","interval")]
colnames(cov_inner)[which(colnames(cov_inner) == "p21003_i2")] <- "age2"
cov_all <- merge(cov_all,cov_inner,by='eid', all = TRUE)

columns_to_standardize <- c('age2', 'Townsend_index', 'BMI2', 'Education','interval','eTIV','FastingTime')

temp <- as.data.frame(scale(cov_all[,..columns_to_standardize]))
cov_all[,columns_to_standardize] = temp[,columns_to_standardize]


# glm regression
list_metafile <- c()
list_targetfile <- c()
list_numberAll <- c()
list_coef <- c()
list_std <- c()
list_tvalues <- c()
list_pvalues <- c()
list_FDR_separate <- c()
FDR_all <- c()

for(i in seq(2,length(meta_index))){

  out1 <- paste0(i-1,"/",length(meta_index)-1)
  print("******************************")
  print(out1)
  print("******************************")
  meta_name <- meta_index[i]
  meta_single <- subset(metadata, select = c(1, i))
  pvalues <- c()

  for(j in seq(2,length(image_index))){
    
    if(j%%10==0){
      out2 <- paste0(j-1,"/",length(image_index)-1)
      print(out2)
    }
    
    img_name <- image_index[j]
    img_single <- subset(imagedata, select = c(1, j))
    list_metafile <- c(list_metafile, meta_name)
    list_targetfile <- c(list_targetfile, img_name)

    useddata <- merge(meta_single, img_single, by=c('eid'))
    useddata <- merge(useddata, cov_all, by=c('eid'))
    useddata <- useddata %>%
      anti_join(ex_df, by = c("eid" = "excluding_edi"))
    

    list_numberAll <- c(list_numberAll, nrow(useddata))

    cleaned_data <- useddata %>% select(meta_name,img_name,"interval","Race_2","Race_3","Race_4","age2","Sex","Townsend_index","BMI2","Education","drinking2_1","drinking2_2","smoking2_1","smoking2_2","Site_11026","Site_11027","eTIV",'FastingTime','Statin')
    colnames(cleaned_data) <- c("meta","IDP","interval","Race_2","Race_3","Race_4","age2","Sex","TDI","BMI2","Education","drinking2_1","drinking2_2","smoking2_1","smoking2_2","Site_11026","Site_11027","eTIV",'FastingTime','Statin')
    
    model <- glm(IDP ~ meta+interval+Race_2+Race_3+Race_4+age2+Sex+TDI+BMI2+Education+smoking2_1+smoking2_2+drinking2_1+drinking2_2+Site_11026+Site_11027+eTIV+FastingTime+Statin, data = cleaned_data)
    summary_model <- summary(model)
    coefficients <- summary_model$coefficients
    
    list_coef <- c(list_coef, coefficients[2,1])
    list_std <- c(list_std, coefficients[2,2])
    list_tvalues <- c(list_tvalues, coefficients[2,3])
    list_pvalues <- c(list_pvalues, coefficients[2,4])
    
    pvalues <- c(pvalues, coefficients[2,4])
    }
  list_FDR_separate <- c(list_FDR_separate,p.adjust(pvalues, method = 'BH'))
}
FDR_all <- p.adjust(list_pvalues,method='BH')

resultframe <- data.frame(meta=list_metafile, img=list_targetfile, numberAll=list_numberAll, 
                          coef=list_coef, std=list_std,tvalue=list_tvalues, 
                          pvalue=list_pvalues, separateFDR=list_FDR_separate,all_FDR=FDR_all)

write.table(resultframe, 'va_glm_meta_IDP.csv',sep=',',row.names = F)
print("Finish!")
