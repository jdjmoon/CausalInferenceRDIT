library(tidyverse)
library(data.table)

cov <- fread('../data/processed_data/processed_covariates_blood_chem_ca.csv')

clinical_trial_dataset <= fread('../data/processed_data/clinical_markers.csv')

protein <- fread('../data/processed_data/processed_blood_chem_ca.csv')
protein <- select(protein, c('eid', 'ADA', 'AGRP', 'ANGPTL3', 'APCS', 'APOA2', 'APOC1', 'APOM', 'CCDC80', 'CCL13', 'CCL27', 'CCN1', 'CES1', 'CREG1', 'DCBLD2', 'ERBB2', 'ERN1', 'FAP', 'FGFBP1', 'GUSB', 'IGSF3', 'ITGB7',
       'ITIH3', 'KLK13', 'LCAT', 'LDLR', 'MSR1', 'NPY', 'PCOLCE', 'PLA2G7', 'PODXL', 'SIGLEC7', 'SMAD5', 'SUMF2', 'TFPI', 'TNFSF10', 'TNFSF11', 'VWA1'))

covs = merge(covs, clinical_trial_dataset)
dt = merge(covs, protein, by = 'eid')

drug = 'C10AA'

cov1 <- paste0("age + gender + as.factor(center) + PC1 + PC2 + bmi + waist_hip_ratio + ", drug)

ct = c('alt', 'ast', 'ggt', 'ck', 'bil', 'alb', 'alp', 'GSI', 'SRM')

cols = colnames(protein)
cols <- cols[-1]

for (j in 1:length(ct)){
    ctype = ct[j]
    print(ctype)
    
    columns = c("index","Estimate","Std..Error", 't.value', 'Pr...t..', 'nrow', 'ctype', 'omics', 'drug') 
    results1 = data.frame(matrix(nrow = 0, ncol = length(columns))) 
    colnames(results1) = columns
    
    for (i in 1:length(cols)) {
        outcome = cols[i]
        print(outcome)
        temp1=na.omit(subset(dt,select=c(outcome,ctype,"age","gender","center","PC1","PC2","bmi","waist_hip_ratio", drug)))
        
        cutoff <- 6
        data <- temp1[[outcome]]
        q1 <- quantile(data, 0.25, na.rm = TRUE)
        q3 <- quantile(data, 0.75, na.rm = TRUE)
        iqr <- q3 - q1
        lower <- q1 - cutoff * iqr
        upper <- q3 + cutoff * iqr
        temp1[[outcome]] <- pmin(pmax(temp1[[outcome]], lower), upper)
        

        result1 = summary(lm(as.formula(parse(text=paste(outcome, "~", ctype, "+", cov1, sep = ""))[[1]]),data=temp1))$coefficients  
        result1 = data.frame(result1)
        result1$nrow = nrow(temp1)
        result1 <- result1 %>% rownames_to_column(var = "index")
        
        result1 = result1[result1$index == ctype | result1$index == drug, ]
        result1$ctype = ctype
        result1$omics = outcome
        result1$drugs = drug
        result1$nlp <- -1 * log(result1$'Pr...t..')
        results1 <- rbind(results1, result1)
    }

    filename <- paste0('../result/statin_prot_', ctype, '.csv')
    fwrite(results1, filename)
}


