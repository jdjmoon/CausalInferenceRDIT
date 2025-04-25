library(tidyverse)
library(data.table)

# Load data
# Blood Biochemistry and Blood Counts = omics_data/bloodchem_firstvisit.csv
# NMR Metabolites = omics_data/ukbnmr_processed_metabolite_firstvisit.csv
# Olink Proteins = omics_data/ukbnmr_processed_metabolite_firstvisit.csv

## Example using the bloodchemistry
omics1 <- fread('../data/omics_data/bloodchem_firstvisit.csv') 
omics1 <- omics1[rowSums(!is.na(omics1)) != 1, ] # Filter rows with only one non-NA value
colnames(omics1)[1] = 'f.eid'
omics2 <- fread('../data/omics_data/bloodchem_secondvisit.csv') 
omics2 <- omics2[rowSums(!is.na(omics2)) != 1, ]
colnames(omics2)[1] = 'f.eid'

# Match by 'f.eid' across both datasets to get common columns
omics1 <-  omics1[omics1$"f.eid" %in%  omics2$"f.eid",]
colnames(omics2) <- gsub("\\.1\\.0$", ".0.0", colnames(omics2))
omics2 <-  omics2[omics2$"f.eid" %in%  omics1$"f.eid",]

# Merge both visits data with suffixes
omics1 <- merge(omics2, omics1, by = "f.eid", suffixes = c(".nmr12", ".nmr11"))

# Apply cutoff for outliers
omicsname <- names(omics1)
omicsname <- omicsname[!omicsname %in% c('f.eid', 'eid', 'visit_index.nmr12', 'visit_index.nmr11')]

for (i in 1:length(omicsname)) {
  outcome <- omicsname[i]
  
  cutoff <- 6
  data <- omics1[[outcome]]
  q1 <- quantile(data, 0.25, na.rm = TRUE)
  q3 <- quantile(data, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - cutoff * iqr
  upper <- q3 + cutoff * iqr
  omics1[[outcome]] <- ifelse(is.na(omics1[[outcome]]), NA,
                                pmin(pmax(omics1[[outcome]], lower), upper))
}

# Compute the difference values between visits
for (col_name in colnames(omics2)) {
  if (col_name != "f.eid") {
    omics1[[col_name]] <- omics1[[paste0(col_name, ".nmr12")]] - omics1[[paste0(col_name, ".nmr11")]]
  }
}

# Select relevant columns and rename
omics <- omics1 %>% select(colnames(omics2))
colnames(omics)[1] = 'f.eid'

# Get list of selected drugs one is Interested in doing analysis
drugs <- c('A10BA', 'A02BC', 'C10AA')

cov <- fread('../data/covariate_related_data/covariates.csv')

# Add blood date to covariate file
blood_date <- fread('../data/covariate_related_data/blood_date.csv')
blood_date <- blood_date %>% select(eid, 'f.3166.0', 'f.3166.1')
blood_date <- na.omit(blood_date)
blood_date$blood_date <- as.Date(substr(blood_date$'f.3166.0', 1, 10))
blood_date$blood_date2 <- as.Date(substr(blood_date$'f.3166.1', 1, 10))
#blood_date$blood_date2 <- as.Date(substr(blood_date$'f.3166.2', 1, 10)) # For the Olink Protein data

path <- '../data/drug_data/ATC4_prescription/'
year <- 3

# Loop through drugs to add 'before' and 'after' prescription status 
for (i in 1:length(drugs)) {
  drug <- drugs[i]
  
  # Before first visit
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
  
  # Before second visit
  name <- paste0(drug, '_af', sep = '')
  file <- fread(paste0(path, drug, '.csv', sep = ''))
  file <- file %>% select(eid, event_dt)
  file$event_dt <- as.Date(file$event_dt)
  file <- merge(blood_date, file, by = 'eid')
  file$tte <- as.integer(file$event_dt - file$blood_date2)
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
fwrite(cov, '../data/processed_data/processed_covariates_blood_chem_la.csv')
fwrite(omics, '../data/processed_data/processed_blood_chem_la.csv')
