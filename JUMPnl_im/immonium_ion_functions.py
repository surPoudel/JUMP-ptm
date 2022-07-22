import pandas as pd
import numpy as np
from all_JUMPnl_functions import ppmCalc,find_prec_match,extract_cols_df


def find_tol(exp_mz_list, intensity_list, imm_ion, noise):
    
    min_tol_val = 10000 # this value will be replaced with the modulus of ppm
    for i,x in enumerate(exp_mz_list):
        if intensity_list[i] > noise:
            tolerance = ppmCalc(imm_ion,x)
            val_tol = abs(tolerance)
            if val_tol < min_tol_val:
                min_tol_val = val_tol
    return min_tol_val

def params_to_dict(parameter): # applied only to NL and immonium ions the parameter should be ion:m/z format separted by ,
    dict1 = {}
    parameter_list = parameter.split(",")
    for imm_ion in parameter_list:
        ion_mz_list = imm_ion.split(":")
        dict1[ion_mz_list[0]] = float(ion_mz_list[1])
    return dict1


def precursor_immonium_loss(df, immonium_ions_dict, tol, noise): # we add nl_list that consists of neutral loss

    # we get intensity controled m/z and intensity and peptide neutral mass
    mz_cols = list(df.columns)
    np_arr = df.to_numpy()
    for imm_ion_sp in immonium_ions_dict.keys():
        imm_ion = immonium_ions_dict[imm_ion_sp.lower()]
        imm_ion_tol_list = []
        for row in np_arr:
            int_list = row[mz_cols.index("int_list_noise_removed")]
            mz_list = row[mz_cols.index("mz_list_noise_removed")]
            tolerance = find_tol(mz_list, int_list, imm_ion, noise)
            match_immion = find_prec_match([imm_ion], mz_list, tol)
            imm_ion_tol_list.append(match_immion)
            
        df[imm_ion_sp+"_lowest_tolerance"] = imm_ion_tol_list


def precursor_immo_count(df,pattern_str="_lowest_tolerance"):
    
    mz_cols = list(df.columns)
    np_arr = df.to_numpy()
    #extract NL_cols
    nl_cols = extract_cols_df(df, pattern_str)
    
    store_nl_all = []
    
    for row in np_arr:
        nl_row_list = []
        for nl in nl_cols:
            nl_score = row[mz_cols.index(nl)]
            if nl_score == 1:
                get_nl = nl.split(pattern_str)[0].upper()
                nl_row_list.append(get_nl)
        
        if len(nl_row_list) == 0:
            outcome_val = "No immonium ion"
        else:
            outcome_val = "/".join(nl_row_list)
        
        store_nl_all.append(outcome_val)
    df["outcome_immonium_ion"] = store_nl_all
    