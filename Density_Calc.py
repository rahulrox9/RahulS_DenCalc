import os
import pandas as pd
import numpy as np


# -------------------------------
# Step 1: Liquid Density Function
# -------------------------------
def calculate_liquid_density(liq, mat):
    params = {
        'SiO2': (60.0855, 26.86, 0.0, -0.000189, 1773),
        'TiO2': (79.88, 28.32, 0.00724, -0.000231, 1773),
        'Al2O3': (101.96, 37.42, 0.00262, -0.000226, 1773),
        'FeO': (71.85, 12.68, 0.00369, -0.000045, 1723),
        'MgO': (40.3, 12.02, 0.00327, 0.000027, 1773),
        'CaO': (56.08, 16.90, 0.00374, 0.000034, 1773),
        'Na2O': (61.98, 29.65, 0.00768, -0.00024, 1773),
        'K2O': (94.2, 47.28, 0.01208, -0.000675, 1773),
        'H2O': (18.02, 22.9, 0.0095, -0.00032, 1273),
    }
    oxides = list(params.keys())
    
    # H2O = MgO * 0.052 
    mat['H2O'] = mat['MgO'] * 0.052 
    liq['H2O'] = liq['MgO'] * 0.052
    
    def compute_density(df, H2O_col):
        df = df.copy()
        df['H2O'] = df[H2O_col]
        df['T_K'] = df['T_degC'] + 273.15
        df['OriginalSum'] = df[oxides].sum(axis=1)
        df[oxides] = df[oxides].div(df['OriginalSum'], axis=0) * 100

        # wt% → mol
        for ox in oxides:
            MW = params[ox][0]
            df[ox + '_mol'] = df[ox] / MW

        # Partial molar volumes
        V_components = []
        for ox in oxides:
            MW, MV, dVdT, dVdP, Tref = params[ox]
            mols = df[ox + '_mol']
            V_partial = MV + dVdT * (df['T_K'] - Tref) + dVdP * (df['P_kbar']*1e3 - 1)
            V_components.append(V_partial * mols)
        df['Vsum'] = sum(V_components)

        # Back to g
        for ox in oxides:
            MW = params[ox][0]
            df[ox + '_g'] = df[ox + '_mol'] * MW
        df['XMW_sum'] = df[oxides].sum(axis=1)
        return df['XMW_sum'] / df['Vsum']

    density_df = pd.DataFrame(index=liq['Sample'].unique())
    density_df['Carrier Melt'] = compute_density(liq, 'H2O').values
    density_df['Interstitial Melt'] = compute_density(mat, 'H2O').values


    return density_df


# -------------------------------
# Step 2: Mineral Density Function
# -------------------------------
def calculate_mineral_density(xl, n_list):
    par_list = ['Ab','An','Or','Fo','Fa','Wo','En','Fs']
    param = {
        'Ab': {'v1': -1.95E-06, 'v2': 0, 'v3': 2.50E-05, 'v4': 6.72E-09, 'ρr': 2.621},
        'An': {'v1': -1.27E-06, 'v2': 3.18E-12, 'v3': 1.09E-05, 'v4': 4.20E-09, 'ρr': 2.765},
        'Or': {'v1': -1.81E-06, 'v2': 5.11E-12, 'v3': 1.51E-05, 'v4': 5.49E-09, 'ρr': 2.571},
        'Fo': {'v1': -7.91E-07, 'v2': 1.35E-12, 'v3': 2.65E-05, 'v4': 8.86E-09, 'ρr': 3.227},
        'Fa': {'v1': -7.30E-07, 'v2': 0, 'v3': 2.65E-05, 'v4': 7.95E-09, 'ρr': 4.402},
        'Wo': {'v1': -1.25E-06, 'v2': 3.11E-12, 'v3': 2.82E-05, 'v4': 0, 'ρr': 2.937},
        'En': {'v1': -7.50E-07, 'v2': 4.48E-13, 'v3': 1.42E-05, 'v4': 3.64E-08, 'ρr': 3.204},
        'Fs': {'v1': -9.90E-07, 'v2': 0, 'v3': 3.18E-05, 'v4': 7.59E-09, 'ρr': 4.002}
    }

    conds = xl.groupby('Sample')[['T_degC','P_kbar']].first().reindex(n_list)
    Tinp = pd.DataFrame(np.tile((conds['T_degC']+273.15).to_numpy().reshape(-1,1), len(par_list)), columns=par_list, index=conds.index)
    Pinp = pd.DataFrame(np.tile((conds['P_kbar']*1000).to_numpy().reshape(-1,1), len(par_list)), columns=par_list, index=conds.index)

    # Molar masses
    MM = { 'FeO': 71.844, 'MgO': 40.3, 'CaO': 56.0774, 'Na2O': 61.979, 'K2O': 94.196}

    # Compute moles of cations
    xl['Ca_apfu'] = (xl['CaO'] / MM['CaO']) * 1 
    xl['Na_apfu'] = (xl['Na2O'] / MM['Na2O']) * 2
    xl['K_apfu'] = (xl['K2O'] / MM['K2O']) * 2   
    xl['Mg_apfu'] = (xl['MgO'] / MM['MgO']) * 1
    xl['Fe_apfu'] = (xl['FeO'] / MM['FeO']) * 1

    # Plagioclase endmembers
    mask_plag = xl['Phase']=='Plagioclase'
    xl.loc[mask_plag, 'An'] = xl.loc[mask_plag, 'Ca_apfu'] / (xl.loc[mask_plag, ['Ca_apfu','Na_apfu','K_apfu']].sum(axis=1))
    xl.loc[mask_plag, 'Ab'] = xl.loc[mask_plag, 'Na_apfu'] / (xl.loc[mask_plag, ['Ca_apfu','Na_apfu','K_apfu']].sum(axis=1))
    xl.loc[mask_plag, 'Or'] = xl.loc[mask_plag, 'K_apfu'] / (xl.loc[mask_plag, ['Ca_apfu','Na_apfu','K_apfu']].sum(axis=1))
    
    # Olivine endmembers
    mask_ol = xl['Phase']=='Olivine'
    xl.loc[mask_ol, 'Fo'] = xl.loc[mask_ol, 'Mg_apfu'] / (xl.loc[mask_ol, ['Mg_apfu','Fe_apfu']].sum(axis=1))
    xl.loc[mask_ol, 'Fa'] = xl.loc[mask_ol, 'Fe_apfu'] / (xl.loc[mask_ol, ['Mg_apfu','Fe_apfu']].sum(axis=1))
    
    # Clinopyroxene endmembers
    mask_cpx = xl['Phase']=='Clinopyroxene'
    xl.loc[mask_cpx, 'Wo'] = xl.loc[mask_cpx, 'Ca_apfu'] / (xl.loc[mask_cpx, ['Ca_apfu','Mg_apfu','Fe_apfu']].sum(axis=1))
    xl.loc[mask_cpx, 'En'] = xl.loc[mask_cpx, 'Mg_apfu'] / (xl.loc[mask_cpx, ['Ca_apfu','Mg_apfu','Fe_apfu']].sum(axis=1))
    xl.loc[mask_cpx, 'Fs'] = xl.loc[mask_cpx, 'Fe_apfu'] / (xl.loc[mask_cpx, ['Ca_apfu','Mg_apfu','Fe_apfu']].sum(axis=1))
    
    min_comp = xl.set_index('Sample')[par_list].groupby('Sample').median()
    den = pd.DataFrame(index=Tinp.index, columns=Tinp.columns)

    for par in par_list:
        p = param[par]
        den[par] = p['ρr'] / (1 + p['v1']*(Pinp[par]-1) + p['v2']*(Pinp[par]-1)**2 + p['v3']*(Tinp[par]-298.15) + p['v4']*(Tinp[par]-298.15)**2)

    result = pd.DataFrame(index=n_list)
    result['Plagioclase'] = den['An']*min_comp['An'] + den['Ab']*min_comp['Ab'] + den['Or']*min_comp['Or']
    result['Olivine'] = den['Fo']*min_comp['Fo'] + den['Fa']*min_comp['Fa']
    result['Clinopyroxene'] = den['Wo']*min_comp['Wo'] + den['En']*min_comp['En'] + den['Fs']*min_comp['Fs']

    return result


# -------------------------------
# Step 3: Bulk Density
# -------------------------------
def compute_bulk_density(density_df, propn):
    return (density_df * propn / 100).sum(axis=1)


# -------------------------------
# Main Workflow
# -------------------------------

data_folder = 'data'
data = pd.read_csv(os.path.join(data_folder,'Average_Phase_Comp.csv'))
PT_df = pd.read_csv(os.path.join(data_folder,'OPAM_HS24_PT.csv'))
propn = pd.read_csv(os.path.join(data_folder,'Modal_Propn.csv'), index_col=0)
propn = propn.rename(columns={'Matrix': 'Interstitial Melt'})

xl = data[data['Sample'] != 'GO19-02.S2']
n_list = xl['Sample'].unique()

PT_df = PT_df.rename(columns={'SAMPLE':'Sample'})[['Sample','T_degC','P_kbar']]
PT_df = PT_df[PT_df['Sample'] != 'GO19-02.S2'].set_index('Sample')

# Prepare matrix and liquid
mat = xl[xl['Phase']=='Matrix'].set_index('Sample').join(PT_df).reset_index()
xl = xl[xl['Phase']!='Matrix'].set_index('Sample').join(PT_df).reset_index()
liq = pd.concat([data[data['Sample']=='GO19-02.S2']]*len(n_list), ignore_index=True)
liq['Sample'] = n_list
liq['Phase'] = 'Liquid'
liq = liq.merge(PT_df.reset_index(), on='Sample', how='left')

# Step 1: Liquid densities
density_liquid = calculate_liquid_density(liq, mat)

# Step 2: Mineral densities
density_mineral = calculate_mineral_density(xl, n_list)

# Combine
Density = pd.concat([density_liquid, density_mineral], axis=1)

# Step 3: Bulk density
Density['Nodule'] = compute_bulk_density(Density[['Carrier Melt','Interstitial Melt','Olivine','Clinopyroxene','Plagioclase']], propn)

os.makedirs('exports', exist_ok=True)
out_file = os.path.join('exports','Density.csv')
Density.to_csv(out_file)
print(f'Density calculations saved to {out_file}')
