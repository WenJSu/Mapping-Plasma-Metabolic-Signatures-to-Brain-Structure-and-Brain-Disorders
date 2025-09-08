library(TwoSampleMR)
library(data.table)

use_data <- read.csv("BONsig_bing_ill.csv")
meta_list <- unique(use_data$meta)
diseases_list <- unique(use_data$diseases)

P_M_MR_hetero <- data.frame()
P_M_MR_pleio <- data.frame()
pro_ukb_IV <- read.csv("cispQTL_IV_UKB.txt",sep='\t')

ins <- format_data(pro_ukb_IV,type = "exposure",header = T,phenotype_col = "Pro_code",snp_col = "rsid",beta_col = "BETA",se_col = "SE",eaf_col = "A1FREQ",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",pval_col = "P")

for (i in 1:length(diseases_list)) {
  output <- paste0(i,"/",length(diseases_list))
  print(output)
  
  ill = diseases_list[i]
  pheno <- as.data.frame(fread(paste0('finngen_R10_',ill),header = T))
  pheno$PHENO <- diseases_list[i]
  
  out <- format_data(pheno,type = "outcome",header = T,snps = ins$SNP,phenotype_col = "PHENO",snp_col = "rsids",beta_col = "beta",se_col = "sebeta",effect_allele_col ="alt",other_allele_col = "ref",eaf_col = "af_alt",pval_col = "mlogp",log_pval = TRUE)
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
  write.csv(results_OR,paste0(diseases_list[i],".csv"),row.names = F)
}

write.csv(P_M_MR_pleio, "all_pleiotropy_tests.csv", row.names = F)
write.csv(P_M_MR_hetero, "all_heterogeneity_tests.csv", row.names = F)
