library(data.table)
library(TwoSampleMR) 
library(ggplot2)
library(tidyverse)
library(ieugwasr)
library(dplyr)
set.seed(2024)

# Read Statin file 
## Use ACE_Paper_GWAS_5.txt if using Ace Inhibitor
file_exp<-'Statin_GWAS_0806.txt' 
ldl_exp_dat_org <- read_exposure_data( 
    filename = file_exp,
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "Freq.A1.1000G.EUR",
    pval_col = "P-value"
)

##################################
# Read Statin Outcome File
## Use OUT_ACEi_SNP_0813.csv if using Ace Inhibitor
file_proteins<-"OUT_Statin_SNP_NMR_0806.csv"
protein_all <- fread(file_proteins, header=TRUE, sep=",")
# Protein_Name<-strsplit(protein_all$FILE,"_")
Protein_Name_List<-rep("",max(protein_all$FILE_Index))
for (i in 1:length(protein_all$FILE_Index)) {
    idx <- protein_all$FILE_Index[i]
    Protein_Name_List[idx] <- protein_all$FILE[i]
}

outcome_dat <- read_outcome_data(
    filename = file_proteins,
    sep = ",",
    snp_col = "ID",
    beta_col = "ES",
    se_col = "SE",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "AF",
    pval_col = "P",
    phenotype_col = "FILE_Index"
)


# Harmonize
dat <- harmonise_data(
    exposure_dat = ldl_exp_dat_org, 
    outcome_dat = outcome_dat
)
dat$protein = Protein_Name_List[dat$outcome]



############################## Statin MR 
# Perform TwoSampleMR using Weighted Median Method
res_wMedian <- mr(dat, method_list = c( "mr_weighted_median"))
res_wMedian$ProteinName<-Protein_Name_List[res_wMedian$outcome]

# Perform Single SNP MR If interested in single snp analysis
res_single <- mr_singlesnp(dat, all_method="mr_weighted_median")
write.csv(res_single, "Statin_SingleSNP.csv", row.names=FALSE)

write.table(res_wMedian, file="Statin_wMedian_NMR_0806.csv", col.names=TRUE, row.names=FALSE, sep=",")
