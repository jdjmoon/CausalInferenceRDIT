# run_model.R

library(tidyverse)
library(data.table)

# Load processed data
cov <- fread('../data/processed_data/processed_covariates_blood_chem_la.csv')
protein <- fread('../data/processed_data/processed_blood_chem_la.csv')

drug_list <- paste(c('A10BA', 'C03AA', 'A02BC', 'B01AA', 'C09AA', 'B01AC', 
                     'H03AA', 'C10AA', 'C07AB', 'C08CA', 'N02BE', 'C09CA', 
                     'J01CA', 'M01AE', 'J01CF', 'M01AB', 'J01FA', 'N02AA', 
                     'R01AD', 'J01AA', 'M02AA', 'N06AA', 'H02AB', 'R03AC'), "_bf", sep = "")
drug_list2 <- paste(c('A10BA', 'C03AA', 'A02BC', 'B01AA', 'C09AA', 'B01AC', 
                      'H03AA', 'C10AA', 'C07AB', 'C08CA', 'N02BE', 'C09CA', 
                      'J01CA', 'M01AE', 'J01CF', 'M01AB', 'J01FA', 'N02AA', 
                      'R01AD', 'J01AA', 'M02AA', 'N06AA', 'H02AB', 'R03AC'), "_af", sep = "")
drug_list3 <- c('A10BA', 'C03AA', 'A02BC', 'B01AA', 'C09AA', 'B01AC', 'H03AA', 'C10AA', 'C07AB', 
                'C08CA', 'N02BE', 'C09CA', 'J01CA', 'M01AE', 'J01CF', 'M01AB', 'J01FA', 'N02AA', 
                'R01AD', 'J01AA', 'M02AA', 'N06AA', 'H02AB', 'R03AC')

columns <- c("index", "Estimate", "Std..Error", 't-value', 'Pr...t..', 'nrow', 'omics', 'drugs')
results1 <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(results1) = columns

cov1 <- "age + gender + as.factor(center) + PC1 + PC2 + bmi + waist_hip_ratio"

for (k in 1:length(drug_list)) {
    print(drug_list3[k])
    # Assign 1 if did not have prescriptions
    # Assign 2 if had prescription continuosuly
    # Assign 3 if started the prescription after the first visit
    # Assign 4 if stopped the prescription after the first visit
    cov$status <- ifelse(cov[[drug_list[k]]] == 0 & cov[[drug_list2[k]]] == 0, 1,
                         ifelse(cov[[drug_list[k]]] == 1 & cov[[drug_list2[k]]] == 1, 2,
                                ifelse(cov[[drug_list[k]]] == 0 & cov[[drug_list2[k]]] == 1, 3, 4)))

    dt <- merge(cov, protein, by = "f.eid", all.x = TRUE)
    cols <- colnames(protein)
    cols <- cols[-1]

    for (i in 1:length(cols)) {
      outcome <- cols[i]
      print(outcome)
        if (drug_list3[k] == 'A10BA') {
            temp1=na.omit(subset(dt,select=c(outcome,drug_list[k],"age","gender","center","PC1","PC2","edu_yrs","tdi","bmi","waist_hip_ratio","current_smk","log_alcohol_weekly_g", 'A02BC_bf', 'C10AA_bf', 'status')))
            cov1="age + gender + as.factor(center) + PC1 + PC2 + as.factor(edu_yrs) + tdi + bmi + waist_hip_ratio + current_smk + log_alcohol_weekly_g + A02BC_bf + C10AA_bf"
        } else if (drug_list3[k] == 'A02BC') {
          temp1=na.omit(subset(dt,select=c(outcome,drug_list[k],"age","gender","center","PC1","PC2","edu_yrs","tdi","bmi","waist_hip_ratio","current_smk","log_alcohol_weekly_g", 'A10BA_bf', 'C10AA_bf', 'status')))
            cov1="age + gender + as.factor(center) + PC1 + PC2 + as.factor(edu_yrs) + tdi + bmi + waist_hip_ratio + current_smk + log_alcohol_weekly_g + A10BA_bf + C10AA_bf"
        } else if (drug_list3[k] == 'C10AA') {
            temp1=na.omit(subset(dt,select=c(outcome,drug_list[k],"age","gender","center","PC1","PC2","edu_yrs","tdi","bmi","waist_hip_ratio","current_smk","log_alcohol_weekly_g", 'A10BA_bf', 'A02BC_bf', 'status')))
            cov1="age + gender + as.factor(center) + PC1 + PC2 + as.factor(edu_yrs) + tdi + bmi + waist_hip_ratio + current_smk + log_alcohol_weekly_g + A10BA_bf + A02BC_bf"
        } else{
            temp1=na.omit(subset(dt,select=c(outcome,drug_list[k],"age","gender","center","PC1","PC2","edu_yrs","tdi","bmi","waist_hip_ratio","current_smk","log_alcohol_weekly_g", 'A10BA_bf', 'A02BC_bf', 'C10AA_bf', 'status')))
            cov1="age + gender + as.factor(center) + PC1 + PC2 + as.factor(edu_yrs) + tdi + bmi + waist_hip_ratio + current_smk + log_alcohol_weekly_g + A10BA_bf + A02BC_bf + C10AA_bf"
        } 
      temp1$status <- as.factor(temp1$status)
      temp1$status <- relevel(temp1$status, ref = 1)

      result1 <- summary(lm(as.formula(paste(outcome, "~ status +", cov1)), data = temp1))$coefficients
      result1 <- data.frame(result1)
      result1$nrow <- nrow(temp1)
      result1 <- result1 %>% rownames_to_column(var = "index")
      result1$omics <- outcome
      result1$drugs <- drug_list3[k]
      results1 <- rbind(results1, result1)
    }
}

# Save the model results
fwrite(results1, '../result/longitudinal_association_model_blood_chem.csv')


