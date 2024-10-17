import numpy as np
import pandas as pd
import math
import os

from rdrobust import rdrobust,rdbwselect,rdplot

from numpy.linalg import LinAlgError

import sys

def rd(drugs):

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

            df = pd.DataFrame(columns=['Coeff', 'Std. Err.', 't-stat.', 'P>|t|', 'CI Lower', 'CI Upper',
            'kernel', 'Omics', 'types','covariates','n','n1','n2', 'bandwidth'])

            for j in dt2_original.columns[1:-68]:
                dt2_original2 = dt2_original.copy()
                y = dt2_original2[j]
                x = dt2_original2.tmp_tte2
            
                covs = dt2_original2[covariates]
                for k in range(len(bw)):
                    for tp in range(len(ker)):
                        try:
                            est1 = rdrobust(y=y, x=x, covs=covs, bwselect = bw[k], kernel = ker[tp])
                            result = est1.coef.join(est1.se, how='outer')
                            result = result.join(est1.t, how='outer')
                            result = result.join(est1.pv, how='outer')
                            result = result.join(est1.ci, how='outer')
                            result['kernel'] = est1.kernel
                            result['Omics'] = j
                            result = result.reset_index()
                            result = result.rename(columns={'index': 'types'})
                            cols = result.columns.tolist()
                            cols = cols[1:] + cols[:1]
                            result = result[cols]
                            result['covariates'] = '+'.join(covariates)
                            result['n'] = np.sum(est1.N)
                            result['n1'] = est1.N[0]
                            result['n2'] = est1.N[1]
                            result['bandwidth'] = bw[k]
                            #df = df.append(result)
                            df = pd.concat([df, result])
                        except:
                            print(f"Error occurred while processing {drug} - {j} with covariate")
                            result = pd.DataFrame(columns=['Coeff', 'Std. Err.', 't-stat.', 'P>|t|', 'CI Lower', 'CI Upper',
                    'kernel', 'Omics', 'types','covariates', 'bandwidth'])
                            result['Omics'] = j
                            df = pd.concat([df, result])
                            not_working.append(drug +'='+j)
            folder = '../result/RDiT_Result/'
            if not os.path.exists(folder):
                # Create a new folder
                os.makedirs(folder)
            csv_name = folder + drug  + '_w_covaraite_B'+'.csv' # Make _M if Metabolite and _P if Protein
            df.to_csv(csv_name)
            
            
            
    print('Finished')

    print(not_working)

if __name__ == "__main__":
    drugs = sys.argv[1:]
    rd(drugs)