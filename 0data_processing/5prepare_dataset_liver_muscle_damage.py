import pandas as pd
import numpy as np
import os

from datetime import timedelta

# create dataframe for each clinical markers
def extract_marker(df2, df3, read2_list, read3_list, marker_name, output_dir="../data/processed_data"):
    df2_marker = df2[df2['read_2'].isin(read2_list)]
    df3_marker = df3[df3['read_3'].isin(read3_list)]
    df_combined = pd.concat([df2_marker, df3_marker]).reset_index(drop=True)
    output_path = os.path.join(output_dir, f"all_{marker_name.lower()}.csv")
    df = df_combined.copy()
    df['value1'] = pd.to_numeric(df['value1'], errors='coerce')
    df['value'] = df['value1']
    
    if 'value2' in df.columns:
        df['value2'] = pd.to_numeric(df['value2'], errors='coerce')
        df['value'] = df['value'].fillna(df['value2'])
    if 'value3' in df.columns:
        df['value3'] = pd.to_numeric(df['value3'], errors='coerce')
        df['value'] = df['value'].fillna(df['value3'])
    df = df[['eid', 'event_dt', 'value']]
    df['value'] = df['value'].replace(0, np.nan)
    df = df.dropna(subset = ['eid', 'event_dt', 'value']).reset_index(drop = True)
    df = df.rename(columns = {'event_dt': 'ct_date'})
    df.set_index('eid').to_csv(output_path)
    return df

# create datafraem to merge the clinical marker to the clcinical dataframe
def merge_lab_within_window(clin_df, statin_df, lab_df, lab_name):
    if lab_name.lower() == 'ck':
        gender_map = clin_df.set_index('eid')['gender']
        lab_df['gender'] = lab_df['eid'].map(gender_map)

        # Compute flags
        lab_df['GSI'] = np.where(
            ((lab_df['gender'] == 1) & (lab_df['value'] >= 180)) |
            ((lab_df['gender'] == 0) & (lab_df['value'] >= 120)),
            1, 0
        )
        lab_df['SRM'] = np.where(
            ((lab_df['gender'] == 1) & (lab_df['value'] >= 180 * 4)) |
            ((lab_df['gender'] == 0) & (lab_df['value'] >= 120 * 4)),
            1, 0
        )
        lab_df = lab_df.drop(['gender'], axis = 1)

    df1 = pd.merge(statin_df, lab_df, on='eid', how='left')
    df1 = df1[(df1['start_date'] <= df1['ct_date']) & (df1['ct_date'] <= df1['end_date'])]

    if lab_name.lower() == 'ck':
        df1a = df1[['eid', 'GSI', 'SRM']]

    df1 = df1.sort_values(['eid', 'ct_date']).groupby('eid').last().reset_index()[['eid', 'value']]
    df2 = pd.merge(clin_df[['eid', 'lower_bound', 'blood_date']], lab_df, on='eid', how='left')
    df2 = df2[~df2['eid'].isin(df1['eid'])]
    df2 = df2[(df2['lower_bound'] <= df2['ct_date']) & (df2['ct_date'] <= df2['blood_date'])]

    if lab_name.lower() == 'ck':
        df2a = df2[['eid', 'GSI', 'SRM']]
        dfat = pd.concat([df1a, df2a])
        dfat1 = dfat.groupby('eid')['GSI'].max().reset_index()
        dfat2 = dfat.groupby('eid')['SRM'].max().reset_index()
        dfat = pd.merge(dfat1, dfat2, on = 'eid')

    df2 = df2.sort_values(['eid', 'ct_date']).groupby('eid').last().reset_index()[['eid', 'value']]

    df_combined = pd.concat([df1, df2])
    df_combined.columns = ['eid', lab_name]

    merged = pd.merge(clin_df, df_combined, on='eid', how='left')
    if lab_name.lower() == 'ck':
        merged = pd.merge(merged, dfat, on = 'eid', how = 'left')

    return merged

# main function to create into single code
def create_clinical_markers():
    ca_all = pd.read_csv('../data/processed_data/processed_covariates_blood_chem_ca.csv', low_memory=False)
    if 'f.eid' in ca_all.columns:
        ca_all = ca_all.rename(columns ={'f.eid': 'eid'})
    all_eids = ca_all[['eid']]
    gp_clinical = pd.read_csv('../data/gp_clinical.csv', low_memory=False)
    df = gp_clinical[gp_clinical['eid'].isin(all_eids['eid'])]
    df2 = df.copy()
    df3 = df.copy()

    dfalp = extract_marker(df2, df3, ['44F..', '44F1.', '44F2.', '44F3.', '44F4.', '44FZ.'], ['44F..', '44F1.', '44F2.', '44F3.', '44F4.', '44FZ.', 'X80DP', 'XE2px', ''], 'alp')
    dfck = extract_marker(df2, df3, ['44H4.', '44HE.', '44HG.'], ['44H4.', '44HE.', '44HG.', 'XE28Y', 'XaES5', 'XaESB'], 'ck')
    dfalt = extract_marker(df2, df3, ['44G..', '44G3.', '44G30', '44G31', '44GA.', '44GB.'], ['44G..', '44G3.', '44G30', '44G31', 'XaIRi', 'XaLJx', 'X80DL', 'X771e', 'X80DL', 'XaLJx', 'X771f'], 'alt')
    dfast = extract_marker(df2, df3, ['44H50', '44H51', '44H52', '44H5.', '44HB.'], ['44H50', '44H51', '44H52', '44H5.', '44HB.', 'X771i', 'X80De'], 'ast')
    dfggt = extract_marker(df2, df3, ['44G4.', '44G40', '44G41', '44G7.', '44G9.'], ['XE28U', 'X80E1', '44G4.', '44G40', '44G41', 'XaES4', 'XaES3', 'XM1D6'], 'ggt')
    dfalb = extract_marker(df2, df3, ['44M4.', '44M40', '44M41', '44MI.'], ['44M4.', '44M40', '44M41', 'X772i', 'XE2eA'], 'alb')
    dfbil = extract_marker(df2, df3, ['44E..', '44E1.', '44E2.', '44E3.', '44E6.', '44E9.', '44EZ.', '44EC.'], ['44E..', '44E1.', '44E2.', '44E3.', '44E6.', '44EZ.', '44EC.', 'XE28O', 'XE2qu', 'Xa972', 'XaERu', 'XaETf'], 'bil')

    statin1 = pd.read_csv('../data/drug_data/ATC4_prescription/C10AA.csv', low_memory=False)
    statin1['event_dt'] = pd.to_datetime(statin1['event_dt'], errors='coerce')
    statin1['start_date'] = statin1['event_dt'] - pd.Timedelta(days=90)
    statin1['end_date'] = statin1['event_dt'] + pd.Timedelta(days=7)
    ca_all1 = ca_all.copy()

    bd = pd.read_csv('../data/covariate_related_data/blood_date.csv')
    bd = bd[['eid', 'f.3166.0']]
    bd.columns = ['eid', 'blood_date']
    bd['blood_date'] = pd.to_datetime(bd['blood_date'])
    bd['lower_bound'] = bd['blood_date'] - pd.to_timedelta(365*3+1,unit='days')
    ca_all1 = pd.merge(ca_all1, bd, on = 'eid')
    bd2 = statin1[['eid','event_dt']]
    bd2['C10AA'] = 1
    bds = pd.merge(bd, bd2, on = 'eid')
    bds['event_dt'] = pd.to_datetime(bds['event_dt'])
    bds['blood_date'] = pd.to_datetime(bds['blood_date'])
    bds = bds[bds['event_dt'] <= bds['blood_date']] 
    bds = bds[bds['lower_bound'] <= bds['event_dt']]
    bds = bds[['eid', 'C10AA']].drop_duplicates().reset_index(drop=True)
    ca_all1 = pd.merge(ca_all1, bds, on='eid', how='left')
    ca_all1['C10AA'] = ca_all1['C10AA'].fillna(0)

    ca_all1 = merge_lab_within_window(ca_all1, statin1, dfalt, 'alt')
    ca_all1 = merge_lab_within_window(ca_all1, statin1, dfast, 'ast')
    ca_all1 = merge_lab_within_window(ca_all1, statin1, dfck, 'ck')
    ca_all1 = merge_lab_within_window(ca_all1, statin1, dfggt, 'ggt')
    ca_all1 = merge_lab_within_window(ca_all1, statin1, dfbil, 'bil')
    ca_all1 = merge_lab_within_window(ca_all1, statin1, dfalb, 'alb')
    ca_all1 = merge_lab_within_window(ca_all1, statin1, dfalp, 'alp')

    ca_all1 = ca_all1[['eid', 'C10AA','alt', 'ast', 'ggt', 'ck', 'bil', 'alb', 'alp', 'GSI', 'SRM']].drop_duplicates().reset_index(drop=True)
    ca_all1.set_index('eid').to_csv('../data/processed_data/clinical_markers.csv')

if __name__ == "__main__":
    create_clinical_markers()