library(tidyverse)
library(data.table)

# Load data
# Blood Biochemistry and Blood Counts = omics_data/bloodchem_firstvisit.csv
# NMR Metabolites = omics_data/ukbnmr_processed_metabolite_firstvisit.csv
# Olink Proteins = omics_data/ukbnmr_processed_metabolite_firstvisit.csv

## Example using the bloodchemistry
omics <- fread('../data/omics_data/bloodchem_firstvisit.csv') 
colnames(omics) <- gsub("-", "_", colnames(omics))
colnames(omics)[1] = 'f.eid'

# Get list of selected drugs one is Interested in doing analysis
drugs <- c('A10BA', 'A02BC', 'C10AA')


cov <- fread('../data/covariate_related_data/covariates.csv')

# Add blood date to covariate file
blood_date <- fread('../data/covariate_related_data/blood_date.csv')
blood_date <- blood_date %>% select(eid, 'f.3166.0')
blood_date <- na.omit(blood_date)
blood_date$blood_date <- as.Date(substr(blood_date$'f.3166.0', 1, 10))

path <- '../data/drug_data/ATC4_prescription/'

# Change number of years if interested in looking at different duration
year <- 3

# Loop through drugs and process prescription data to add if individual had prescription within the duration
for (i in 1:length(drugs)) {
  drug <- drugs[i]
  name <- paste0(drug, '_bf', sep = '')
  file <- fread(paste0(path, drug, '.csv', sep = ''))
  file <- file %>% select(eid, event_dt)
  file$event_dt <- as.Date(file$event_dt)
  file <- merge(blood_date, file, by = 'eid')
  file$tte <- as.integer(file$event_dt - file$blood_date)
  file[[name]] <- as.integer((file$tte > (-365 * year - 1)) & (file$tte < 0))

  file <- file %>%
    group_by(eid) %>%
    summarise(max_value = max(!!sym(name), na.rm = TRUE))
  file[[name]] <- file$max_value
  file <- file %>% select(eid, all_of(name))
  colnames(file)[1] = 'f.eid'

  cov <- merge(cov, file, all.x = TRUE, by = 'f.eid')
  cov <- cov %>% mutate(!!name := replace_na(!!sym(name), 0))
}

# Save the processed data for later use in model running
fwrite(cov, '../data/processed_data/processed_covariates_blood_chem_ca.csv')
fwrite(omics, '../data/processed_data/processed_blood_chem_ca.csv')

