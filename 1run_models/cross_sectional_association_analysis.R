library(tidyverse)
library(data.table)

# Load processed data
cov <- fread('../data/processed_data/processed_covariates_blood_chem_ca.csv')
protein <- fread('../data/processed_data/processed_blood_chem_ca.csv')

drug_list <- c('A10BA_bf', 'C03AA_bf', 'A02BC_bf', 'B01AA_bf', 'C09AA_bf', 'B01AC_bf', 
               'H03AA_bf', 'C10AA_bf', 'C07AB_bf', 'C08CA_bf', 'N02BE_bf', 'C09CA_bf', 
               'J01CA_bf', 'M01AE_bf', 'J01CF_bf', 'M01AB_bf', 'J01FA_bf', 'N02AA_bf',
               'R01AD_bf', 'J01AA_bf', 'M02AA_bf', 'N06AA_bf', 'H02AB_bf', 'R03AC_bf')

cov1 <- "age + gender + as.factor(center) + PC1 + PC2 + bmi + waist_hip_ratio"
columns <- c("index", "Estimate", "Std..Error", 't-value', 'Pr...t..', 'nrow', 'drugs', 'omics') 
results1 <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(results1) = columns

dt <- merge(cov, protein, by = "f.eid", all.x = TRUE)
cols <- colnames(protein)[-1]

for (j in 1:length(drug_list)) { 
    drug <- drug_list[j]
    print(drug)
    for (i in 1:length(cols)) {
        outcome <- cols[i]
        print(outcome)
        if (drug == 'A10BA_bf') {
          temp1 <- na.omit(subset(dt, select=c(outcome,drug,"age","gender","center","PC1","PC2","edu_yrs","tdi","bmi","waist_hip_ratio","current_smk","log_alcohol_weekly_g", 'A02BC_bf', 'C10AA_bf')))
        } else if (drug == 'A02BC_bf') {
          temp1 <- na.omit(subset(dt, select=c(outcome,drug,"age","gender","center","PC1","PC2","edu_yrs","tdi","bmi","waist_hip_ratio","current_smk","log_alcohol_weekly_g", 'A10BA_bf', 'C10AA_bf')))
        } else if (drug == 'C10AA_bf') {
            temp1 <- na.omit(subset(dt, select=c(outcome,drug,"age","gender","center","PC1","PC2","edu_yrs","tdi","bmi","waist_hip_ratio","current_smk","log_alcohol_weekly_g", 'A10BA_bf', 'A02BC_bf')))
        } else {
            temp1 <- na.omit(subset(dt, select=c(outcome,drug,"age","gender","center","PC1","PC2","edu_yrs","tdi","bmi","waist_hip_ratio","current_smk","log_alcohol_weekly_g", 'A10BA_bf', 'A02BC_bf', 'C10AA_bf')))
        }

        # Outlier handling
        cutoff <- 6
        data <- temp1[[outcome]]
        q1 <- quantile(data, 0.25, na.rm = TRUE)
        q3 <- quantile(data, 0.75, na.rm = TRUE)
        iqr <- q3 - q1
        lower <- q1 - cutoff * iqr
        upper <- q3 + cutoff * iqr
        temp1[[outcome]] <- pmin(pmax(temp1[[outcome]], lower), upper)

        # Run linear model
        result1 <- summary(lm(as.formula(parse(text=paste(outcome, "~", drug, "+", cov1, sep = ""))[[1]]), data=temp1))$coefficients
        result1 <- data.frame(result1)
        result1$nrow <- nrow(temp1)
        result1 <- result1 %>% rownames_to_column(var = "index")
        result1$drugs <- substr(drug, 1, nchar(drug)-3)
        result1$omics <- outcome
        result1 <- subset(result1, index == drug)
        results1 <- rbind(results1, result1)
    }
}

# Save the model results
fwrite(results1, '../result/cross_sectional_association_model_blood_chem.csv')