#Proteome-Metabolome cis-MR#
library(TwoSampleMR)
library(data.table)

refer_list <- read.csv("s3_del_reverse.csv")
list <- unique(refer_list$meta)

P_M_MR_hetero <- data.frame()
P_M_MR_pleio <- data.frame()
pro_ukb_IV <- read.csv("cispQTL_IV_UKB.txt",sep='\t')
ins <- format_data(pro_ukb_IV,type = "exposure",header = T,phenotype_col = "Pro_code",snp_col = "rsid",beta_col = "BETA",se_col = "SE",eaf_col = "A1FREQ",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",pval_col = "P")

for (i in 1:length(list)) {
  output <- paste0(i,"/",length(list))
  print(output)
  pheno <- data.frame()
  for (j in 1:22) {
    dt <- as.data.frame(fread(paste0(j,".",list[i],".glm.linear"),header = T))
    pheno <- rbind(pheno,dt)
  }
  pheno$PHENO <- list[i]
  out <- format_data(pheno,type = "outcome",header = T,snps = ins$SNP,phenotype_col = "PHENO",snp_col = "ID",beta_col = "BETA",se_col = "SE",effect_allele_col ="A1",other_allele_col = "AX",eaf_col = "A1_FREQ",pval_col = "P")
  harmo <- harmonise_data(exposure_dat = ins,outcome_dat = out)
  mr_res <- mr(harmo,method_list = c("mr_wald_ratio","mr_ivw"))
  mr_pleio <- mr_pleiotropy_test(harmo)
  mr_hetero <- mr_heterogeneity(harmo)
  results_OR <- generate_odds_ratios(mr_res)
  results_OR <- results_OR[,-c(1,2)]
  mr_pleio <- mr_pleio[,-c(1,2)]
  mr_hetero <- mr_hetero[,-c(1,2)]
  P_M_MR_pleio <- rbind(P_M_MR_pleio,mr_pleio)
  P_M_MR_hetero <- rbind(P_M_MR_hetero,mr_hetero)
  write.csv(results_OR,paste0(list[i],".csv"),row.names = F)
}

write.csv(P_M_MR_pleio, "all_pleiotropy_tests.csv", row.names = F)
write.csv(P_M_MR_hetero, "all_heterogeneity_tests.csv", row.names = F)
