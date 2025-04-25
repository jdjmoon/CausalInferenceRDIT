import pandas as pd
import numpy as np
import os

# Get list of selected drugs one is Interested in doing analysis
drugs = ['A10BA', 'A02BC', 'C10AA']

cov = pd.read_csv('../data/covariate_related_data/covariates.csv', low_memory = False)
cov = cov.rename(columns = {'f.eid': 'eid'})

pres = pd.read_csv('../data/drug_data/ATC4_prescription/complete_gp_prescription_with_ATC4.csv', low_memory = False)
pres['event_dt'] = pd.to_datetime(pres['event_dt']).dt.date

# Add blood date to covariate file
blood_date = pd.read_csv('../data/covariate_related_data/blood_date.csv')
blood_date = blood_date[['eid', 'f.3166.0']].dropna().reset_index(drop = True)
blood_date.columns = ['eid', 'blood_date']
blood_date['blood_date'] = pd.to_datetime(blood_date['blood_date']).dt.date

pres2 = pres.copy()
cov2 = cov.copy()
blood_date2 = blood_date.copy()


# Load data
# Blood Biochemistry and Blood Counts = omics_data/bloodchem_firstvisit.csv
# NMR Metabolites = omics_data/ukbnmr_processed_metabolite_firstvisit.csv
# Olink Proteins = omics_data/ukbnmr_processed_metabolite_firstvisit.csv

nmr = pd.read_csv('../data/omics_data/bloodchem_firstvisit.csv')
nmr = nmr.rename(columns = {'f.eid':'eid'})

pres = pres2.copy()
cov = cov2.copy()
blood_date = blood_date2.copy()

cov = cov[cov['eid'].isin(nmr['eid'])].reset_index(drop = True)
nmr = nmr[nmr['eid'].isin(cov['eid'])].reset_index(drop = True)
pres = pres[pres['eid'].isin(nmr['eid'])].reset_index(drop = True)
blood_date = blood_date[blood_date['eid'].isin(nmr['eid'])].reset_index(drop = True)

pres = pd.merge(pres, blood_date, on ='eid')
pres['event_dt'] = pd.to_datetime(pres['event_dt']).dt.date
pres['blood_date'] = pd.to_datetime(pres['blood_date']).dt.date

# Calculate Number of days 
tte2 = (pres['event_dt'] - pres['blood_date'])
tte = []
for i in tte2:
    tte.append(i.days)
pres['tte'] = tte    


# Get first date of prescription for each ATC Level 4 Code
pres['ATC4'] = pres['ATC'].str[:5]
pres = pres[['eid', 'event_dt', 'blood_date', 'ATC', 'ATC4', 'mapping_term', 'tte']]
pres = pres.sort_values(['eid', 'event_dt', 'ATC4']).reset_index(drop = True)
pres = pres.drop_duplicates(subset =['eid', 'ATC4']).reset_index(drop = True)

# Only extract selected drugs
pres = pres[pres['ATC4'].isin(drugs)].reset_index(drop = True)

pres = pd.merge(pres, cov, on='eid')
pres_org = pres

years = [3]

for year in years:
    print(year)
    # Extract the data to be within 3 years from the blood date
    duration=365*year + 1
    duration_y = '_' + str(year) + 'y'
    pres=pres_org
    pres1 = pres[pres['tte'] > 0].reset_index(drop = True)
    pres = pres1
    
    add = pres[pres['ATC4']=='A10BA'].reset_index(drop = True)
    add = add[['eid']]
    add['A10BA'] = 1
    pres = pd.merge(pres, add, on = 'eid', how = 'left')
    pres['A10BA'] = pres['A10BA'].fillna(0)

    add = pres[pres['ATC4']=='C10AA'].reset_index(drop = True)
    add = add[['eid']]
    add['C10AA'] = 1
    pres = pd.merge(pres, add, on = 'eid', how = 'left')
    pres['C10AA'] = pres['C10AA'].fillna(0)

    add = pres[pres['ATC4']=='A02BC'].reset_index(drop = True)
    add = add[['eid']]
    add['A02BC'] = 1
    pres = pd.merge(pres, add, on = 'eid', how = 'left')
    pres['A02BC'] = pres['A02BC'].fillna(0)
    
    for i in range(len(drugs)):
        dt = pres[pres['ATC4']==drugs[i]].reset_index(drop = True)
        name = drugs[i] + duration_y
        dt2 = pd.merge(nmr, dt, on = 'eid')
        if len(dt2) > 1:
            for j in range(1, len(nmr.columns)):
                cutoff = 6
                data = dt2[nmr.columns[j]]
                q3, q1 = np.percentile(data.dropna(), [75 ,25], axis = 0)
                iqr = q3 - q1
                lower = q1 - cutoff*iqr
                upper = q3 + cutoff*iqr
                bools = dt2[nmr.columns[j]].isna()
                dt2[nmr.columns[j]].where(dt2[nmr.columns[j]] <= upper, upper, inplace = True)
                dt2[nmr.columns[j]].where(dt2[nmr.columns[j]] >= lower, lower, inplace = True)
                dt2.loc[bools, nmr.columns[j]] = np.nan
            dt2[name] = np.where(dt2['tte'] < 0, 1, 0)
            dt2['tmp_tte'] = - dt2['tte']
            dt2 = dt2.dropna(subset = ['age', 'gender', 'center', 'PC1', 'PC2', 'bmi',"edu_yrs", "tdi", "waist_hip_ratio", "current_smk", "log_alcohol_weekly_g", 'A10BA', 'C10AA', 'A02BC']).reset_index(drop = True)
            dt2.to_csv('../data/processed_data/DATA_B/new_' + drugs[i] +duration_y+'.csv')
            #table_list.append(dt2)
        else:
            print(name)