#### import statsmodels.formula.api as smf
import numpy as np
import pandas as pd
import math
import os
import seaborn as sns
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

from rdrobust import rdrobust,rdbwselect,rdplot
from rdpower import rdpower, rdsampsi, rdmde
from numpy.linalg import LinAlgError

import sys

def rdpow(drugs):

    cov = pd.read_csv('../data/covariate_related_data/covariates.csv', low_memory = False)
    cov = cov.rename(columns = {'f.eid': 'eid'})
    
    years = [3]
    not_working = []

    covariates = ['age', 'gender', 'center', 'PC1', 'PC2', 'bmi', 'waist_hip_ratio']
    # If testing for different bandwidth, modify bw
    bw = ['mserd']
    # If testing for different kernel,modify ker
    ker = ['uniform']

    for year in years:
        print(year)
        duration_y = '_'+ str(year) + 'y'
        for i in drugs:
            print(i)
            drug = i
            csv_file = 'DATA_B/' + drug  + duration_y+'.csv' # Make _M if Metabolite and _P if Protein


            dt2_original = pd.read_csv(csv_file)
            if 'Unnamed: 0' in dt2_original.columns:
                dt2_original = dt2_original.drop(['Unnamed: 0'], axis = 1)

            drug = drug + duration_y
            dt2_original['tmp_tte2'] = dt2_original['tmp_tte']

            df = pd.DataFrame(columns=['power_rbc', 'se_rbc', 'sampsi_r', 'sampsi_l', 'samph_r', 'samph_l', 'N_r', 'N_l', 'Nh_l', 'Nh_r', 'tau', 'bias_r', 'bias_l', 'Vr_rb', 'Vl_rb', 'alpha', 'Omics', 'drug'])

            for j in dt2_original.columns[1:-68]:
                dt2_original2 = dt2_original.copy()
                Z = dt2_original2[[j, 'tmp_tte2']]
            
                covs = dt2_original2[covariates]
                for k in range(len(bw)):
                    for tp in range(len(ker)):
                        
                        try:
                            aux = rdpower(data=Z, covs = covs, cutoff = 0, kernel = 'uniform', plot = True)
                            
                            tau = aux['tau']
                            result = pd.DataFrame(aux)
                            result['Omics'] = j
                            result['drug'] = i
                            df = pd.concat([df, result])
                            
                            for ts in range(0, 5):
                                aux = rdpower(data=Z, covs=covs, cutoff=0, kernel='uniform', plot=True, tau=0.2 * ts * tau) 
                                result = pd.DataFrame(aux)
                                df = pd.concat([df, result])
                           
                        except:
                            print(f"Error occurred while processing {drug} - {j} with covariate")
                            result = pd.DataFrame(columns=['power_rbc', 'se_rbc', 'sampsi_r', 'sampsi_l', 'samph_r', 'samph_l', 'N_r', 'N_l', 'Nh_l', 'Nh_r', 'tau', 'bias_r', 'bias_l', 'Vr_rb', 'Vl_rb', 'alpha', 'Omics', 'drug'])
                            result['Omics'] = j
                            result['drug'] = i
                            df = pd.concat([df, result])
                            not_working.append(drug +'='+j)
            folder = '../result/RDiT_PowerResult/'
            if not os.path.exists(folder):
                # Create a new folder
                os.makedirs(folder)
            csv_name = folder + drug  + '_w_covaraite_B'+'.csv' # Make _M if Metabolite and _P if Protein
            df.to_csv(csv_name)
            
            
            
    print('Finished')

    print(not_working)

if __name__ == "__main__":
    drugs = sys.argv[1:]
    rdpow(drugs)