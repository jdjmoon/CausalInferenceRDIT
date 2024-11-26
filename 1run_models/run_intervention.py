import numpy as np
import pandas as pd
import math
import os
import seaborn as sns
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

from rdrobust import rdrobust,rdbwselect,rdplot

from numpy.linalg import LinAlgError

import sys

cov = pd.read_csv('../data/processed_data/processed_covariates_blood_chem_ca.csv', low_memory = False)

covariates = ['age', 'gender', 'center', 'PC1', 'PC2', 'bmi', 'waist_hip_ratio']

bw = ['mserd']
ker = ['uniform']

year = 3

# Change Csv file, drug, and omics id based on your need
csv_file = '../data/processed_data/DATA_B/new_A10BA_3y.csv' 
drug = 'A10BA'
omics = 'f.30750.0.0'

dt2_original = pd.read_csv(csv_file)
if 'Unnamed: 0' in dt2_original.columns:
    dt2_original = dt2_original.drop(['Unnamed: 0'], axis = 1)

dt2_original['tmp_tte2'] = dt2_original['tmp_tte']


dt2_original2 = dt2_original.copy()

df = pd.DataFrame(columns=['Coeff', 'Std. Err.', 't-stat.', 'P>|t|', 'CI Lower', 'CI Upper', 'kernel', 'Omics', 'types','covariates','n','n1','n2', 'bandwidth'])

x = dt2_original2.tmp_tte2

np.random.seed(42)

combined_range = np.arange(x.min(), -365*3-1-1)
np.random.shuffle(combined_range)
#randomc = np.random.choice(combined_range)

used_c = set()
successful_runs = 0

while successful_runs < 10 or len(used_c) < len(combined_range):
    randomc = np.random.choice(combined_range)

    if randomc in used_c:
        continue

    try:
        #est1 = rdrobust(y=y, x=x, covs=covs, bwselect=bw[0], kernel=ker[0], c=randomc)
        dt2_original21 = dt2_original2[dt2_original2['tmp_tte'] < randomc + 365*3+1]
        dt2_original21 = dt2_original21[dt2_original21['tmp_tte'] > randomc - 365*3-1].reset_index(drop = True)
        y = dt2_original21[omics]
        x = dt2_original21.tmp_tte2
        covs = dt2_original21[covariates]
        est1 = rdrobust(y=y, x=x, covs=covs, bwselect = bw[0], kernel = ker[0], c= randomc)

        result = est1.coef.join(est1.se, how='outer')
        result = result.join(est1.t, how='outer')
        result = result.join(est1.pv, how='outer')
        result = result.join(est1.ci, how='outer')
        result['kernel'] = est1.kernel
        result['drug'] = drug
        result['Omics'] = omics
        result = result.reset_index()
        result = result.rename(columns={'index': 'types'})

        result['covariates'] = '+'.join(covariates)
        result['n'] = np.sum(est1.N)
        result['n1'] = est1.N[0]
        result['n2'] = est1.N[1]
        result['bandwidth'] = bw[0]
        result['c'] = randomc


        df = pd.concat([df, result])
        used_c.add(randomc)
        successful_runs += 1

    except Exception as e:
        print(f"rdrobust failed for c={randomc}: {e}")
        used_c.add(randomc)
        continue

print(f"Successfully executed rdrobust {successful_runs} times.")


df.set_index('Coeff').to_csv('../result/intervention_'+ drug+'_'+ omics+'.csv')