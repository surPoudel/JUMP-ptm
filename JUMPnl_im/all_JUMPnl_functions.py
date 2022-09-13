import sys
import pandas as pd
import glob
import os
import re
from collections import Counter
import configparser
config = configparser.ConfigParser()

import pyteomics as pyteomics
from pyteomics import mzxml, auxiliary, pepxml, mass, ms2
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
from os.path import dirname
import itertools   

def cpFile(filename, folder):
    cmd1 = "cp {} {}".format(filename,folder)
    os.system(cmd1)

def make_spectrum(row):
    out = row.Outfile
    outfile_jump = os.path.basename(out)
    out_split = outfile_jump.split(".")
    spectrum = out_split[0]+"."+out_split[1]+"."+out_split[-2]
    return spectrum


#check function
def check(filename, fraction):
    with open(filename, "w") as f:
        f.write("{} is analyzed for neutral loss or immonium ion".format(fraction))


#jump_l_functions

#this function extracts dynamic and static modifications using pepxml files and stores them as the dictionary
def getDynStatModsInfoPepXml(pepxml):
    f = open(pepxml,"r") #read the file
    line = f.readline()
    var_AA_mass = {} 
    var_AA_symbol = {} #symbol mass dictionary
    stat_AA_mass = {}
    while "<spectrum_query" not in line: #end reading the file if this is sen
    # look for aminocaid modification keyword so that the modification infomration can be parsed
        if "<aminoacid_modification" in line: #found modification information
            if "symbol=" in line.strip(): #only dynamic modification has symbols as static are fixed

                #check the patter
                pattern = '<aminoacid_modification aminoacid="(\w)" massdiff="([-]?\d+\.\d+)" mass="\d+\.\d+" variable="(\w)" symbol="(.+?)"/>'

                #match for the patter in the line and store them as the tuples there are 4 matches () this is changed to list with list function
                #print ("pattern = ",pattern)
                #print ("   ",line.strip())
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3,4))

                modAA = mods_Description[0] #first element in list is amino acid 
                varMass = float(mods_Description[1]) #second element in list is mass and it is converted to float
                variable = mods_Description[2] #third element in the list is determination of variable of static modification "Y" is variable and "N" is static
                symbol = mods_Description[3] #Symbol. this is used in the dictionary
                valAddKey(var_AA_mass,  varMass, modAA)
    #             valAddKey(var_AA_symbol, symbol, varMass)
                var_AA_symbol[symbol] = varMass #symbol as key and mass as values in the dictionary
                line = f.readline()
            else:
                #this is for static modification so there is no symbol hence we have only three values in the list
                pattern = '<aminoacid_modification aminoacid="(\w)" massdiff="([-]?\d+\.\d+)" mass="\d+\.\d+" variable="(\w)"/>'
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3))
                modAA = mods_Description[0]
                varMass = float(mods_Description[1])
                variable = mods_Description[2]
                stat_AA_mass[modAA] = varMass
    #             valAddKey(stat_AA_mass, modAA, varMass)
                line = f.readline()

        elif "<terminal_modification terminus" in line: #This is for terminal modification such as N-term or C-term
            if "symbol=" in line.strip():
                pattern = '<terminal_modification terminus="(\w)" massdiff="([-]?\d+\.\d+)" mass="\d+\.\d+" variable="(\w)".+symbol="(.+?)"/>'
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3,4))

                modAA = mods_Description[0].lower()
                varMass = float(mods_Description[1])
                variable = mods_Description[2]
                symbol = mods_Description[3]
                valAddKey(var_AA_mass,  varMass, modAA)
    #             valAddKey(var_AA_symbol, symbol, varMass)
                var_AA_symbol[symbol] = varMass
                line = f.readline()
            else:
                pattern = '<terminal_modification terminus="(\w)" massdiff="([-]?\d+\.\d+)" mass="\d+\.\d+" variable="(\w)".+/>'
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3))
                modAA = mods_Description[0].lower()
                varMass = float(mods_Description[1])
                variable = mods_Description[2]
                stat_AA_mass[modAA] = varMass
    #             valAddKey(stat_AA_mass, modAA, varMass)
                line = f.readline()
        else:
            line = f.readline()
    return var_AA_mass,var_AA_symbol, stat_AA_mass

def valAddKey(dict1, key, val):
    if key not in dict1.keys():
    	dict1[key] = [val]
    else:
    	dict1[key].append(val)
    return dict1

# merge all pkl files to make a dataframe

def merge_pkl(pkl_list):
    newDF = pd.read_pickle(pkl_list[0])
    
    for val in range(1, len(pkl_list)):
        newDF_1 = pd.read_pickle(pkl_list[val])
        frames = [newDF, newDF_1]
        newDF = pd.concat(frames)

    return newDF

def mkdir(dir1):
    cmd = "mkdir "+dir1
    try:
        os.system(cmd)
    except:
        "Dynamic Modification Directory exits!!"

# reading jump -f files 

def return_skiprows(file, delimiter, peptide):
    with open(file, "r") as f:
        skiprows = 0
        for line in f:
            if peptide+delimiter in line:
                break
            else:
                skiprows+=1
    return skiprows


def count_rawfile_from_idtxt(idtxt, stage):
   psms_all=pd.read_csv(idtxt, delimiter=";",skiprows=return_skiprows(idtxt, ";", "Peptide"), low_memory=False)
   all_cols = list(psms_all.columns)
   if "ptm_stage" in all_cols:
        psms_all = psms_all.loc[psms_all.ptm_stage == stage]
        if stage == "Stage_6":
            psms_all = psms_all.loc[psms_all.Peptide.str.contains("R@")]

   pho_data_NR = psms_all.drop_duplicates(subset=["Outfile"])
   pho_data_NR["spectrum"] = pho_data_NR.apply(make_spectrum, axis=1)
   pho_data_NR[["exp","scan","charge"]] = pho_data_NR["spectrum"].str.split(".",expand=True)
   final_files = set(pho_data_NR.exp)
       
   return final_files



###        
def read_idtxt(idtxt, ms2File, stage, sta_AA, jump_mod_dict):
    check_variable = os.path.basename(dirname(ms2File))
    psms_all=pd.read_csv(idtxt, delimiter=";",skiprows=return_skiprows(idtxt, ";", "Peptide"), low_memory=False)


    if "Peptides_original" in list(psms_all.columns):
          psms_all["Peptides"] = psms_all.Peptides_original       

    print (".... Reading complete")
    # #check fraction p12
    # psms_all = psms_all.loc[psms_all.Outfile.str.contains("p12")]

    if "lscore" in idtxt:
          psms_all.rename(columns = {"JUMPl_site":"Peptides"}, inplace=True)

    if "Peptides" not in list(psms_all.columns):
          psms_all["Peptides"] = psms_all["Peptide"]

    all_cols = list(psms_all.columns)

    # print (all_cols)
    # print (psms_all.head())

    if "ptm_stage" in all_cols:
          psms_all = psms_all.loc[psms_all.ptm_stage == stage]

    #c = [x for x in a if x not in b]
    remove_cols = ["spectrum","exp","scan","charge","expScan","plain_peptide","modifications","combined"]
    core_cols = [x for x in all_cols if x not in remove_cols]

    psms_all = psms_all[core_cols]
    # print (psms_all.columns)


    pho_data_NR = psms_all.drop_duplicates(subset=["Outfile"])
    pho_data_NR["spectrum"] = pho_data_NR.apply(make_spectrum, axis=1)
    pho_data_NR[["exp","scan","charge"]] = pho_data_NR["spectrum"].str.split(".",expand=True)
    pho_data_NR = pho_data_NR.loc[pho_data_NR["exp"] == check_variable]
    print ("the length of the dataframe is {}".format(pho_data_NR.shape[0]))
    pho_data_NR["expScan"] = pho_data_NR.exp+"."+pho_data_NR.scan.astype("str")
    # print (set(psms_all["ptm_stage"]))

    pho_data_NR[["plain_peptide","modifications"]] = pho_data_NR.apply(computeModifications, sta_AA=sta_AA,jump_mod_dict=jump_mod_dict, axis=1)
    pho_data_NR["combined"] = pho_data_NR["spectrum"]+"\t"+pho_data_NR["plain_peptide"]+"\t"+pho_data_NR["modifications"]
    # pho_data_NR["exp"] = pho_data_NR.Outfile.apply(lambda x: os.path.basename(dirname(dirname(x))))



    required_cols = ['Peptides', 'Outfile','XCorr','spectrum', 'plain_peptide', 'modifications', 'combined','exp','scan','charge','expScan']


    pho_data_NR = pho_data_NR[required_cols]
    # pho_data_NR["jump_ms2"] = pho_data_NR.Outfile.apply(lambda x: dirname(x).split(".")[0]+".ms2")

    # pho_data_NR["exp"] = pho_data_NR.jump_ms2.apply(lambda x: os.path.basename(x))
    pho_data_NR["absolute_presearch_ms2_path"] = ms2File
       
    return pho_data_NR


#This fucntion computes modifications with all the static and dynamic represenatation
def computeModifications(row, jump_mod_dict, sta_AA, peptide = "Peptides"):
    mod_peptide = row[peptide] #Peptide sequence predicted by the software
    peptideNoFlank = mod_peptide.split(".")[1] #Flanking peptide
    pepNoFlankList = list(peptideNoFlank) #flanking peptide is converted to the list, so the symbols are also member of list along with amino acids

    plain_peptide_list =[] #iniiates plain peptide formation
    for aa in pepNoFlankList:
        if aa in pyteomics.parser.std_amino_acids: #this looks at the standard amino acids in pyteomics and if the value in list is not amino acid for example * or other symbols, then it will discard those and just adds amino acdis
            plain_peptide_list.append(aa) 
    plain_peptide = "".join(plain_peptide_list) #creates a new string or amino acids
    
    dynAA_count = 0 #counts the number of dynamic modifcaiton in the peptide sequence
    for key in jump_mod_dict.keys(): 
        if key in mod_peptide:
            dynAA_count +=1 #updates dynamic modification if it finds the symbol in modified peptide

    mods = [] #this initiates the peptide modificaion
    if "n" in peptideNoFlank: #looks for dynamic n-term
        dynAA_count-=1 #this updates dynAA_count to 1 value less if n-term is tehre because n term is already dealt here so we can reduce that by 1
        mods.append("1_V_"+str(jump_mod_dict["n"])+"_n") #comet modificaiotn type generation
        peptideNoFlank = peptideNoFlank[1:] #if n-term is found this will now remove n-term from plain peptide for consideration because n-term is already considered
    #find index of other dynamic modification symbols

    mod_symbol_index = [] #initiates a list with index of modified symbols
    for key in jump_mod_dict.keys(): 
        for aa_index, aa in enumerate(peptideNoFlank): #looks for index and amino acid in pepptide Nofranl (it has symbols)
            if aa == key: #if symbol is found
                mod_symbol_index.append(aa_index) #the index is appended
    no_remaining_dyn_mods = len(mod_symbol_index) #lenght of total dynamic modificaiton except n-term

    mod_symbol_index.sort()  #sorts the index in ascending order so that we get correct modifications sites. This is very important when we have more than one dynamic modiification
    # print ("mod symbol _index =", mod_symbol_index)
    iterate = 1
    for mod_no in range(0,no_remaining_dyn_mods):
        aa = peptideNoFlank[mod_symbol_index[mod_no]-1]
        position = mod_symbol_index[mod_no]-iterate+1 #position is always greater than index
        mod_mass = jump_mod_dict[peptideNoFlank[mod_symbol_index[mod_no]]]

        iterate+=1
        new_mod = str(position)+"_V_"+str(mod_mass)
        mods.append(new_mod)
    #static modification

    for sta_key in sta_AA.keys():
        if sta_key == "n":
            sta_mod = "1_S_"+str(sta_AA["n"])+"_n"
            mods.append(sta_mod)
        else:
            plain_pep_list = list(plain_peptide)
            for index, aa in enumerate(plain_pep_list):
                if aa == sta_key:
                    sta_mod = str(index+1)+"_S_"+str(sta_AA[aa])
                    mods.append(sta_mod)
    modifications = ",".join(mods)  
    # print (modifications)            
    return pd.Series([plain_peptide,modifications])

def cleaning_spectra(df, intensity_cut=0.33, fixed_noise=5000):
    mz_cols = list(df.columns)
    np_arr = df.to_numpy()
    new_mz_list_noise_removal = []
    new_int_list_noise_removal = []
    
    for row in np_arr:

        exp_mz_list = row[mz_cols.index("exp_mz_list")]
        intensity_list = row[mz_cols.index("intensity_list")]
        int_list,mz_list = noise_removal(exp_mz_list, intensity_list, intensity_cut, fixed_noise)
        
        new_mz_list_noise_removal.append(mz_list)
        new_int_list_noise_removal.append(int_list)
        
    df["mz_list_noise_removed"] = new_mz_list_noise_removal
    df["int_list_noise_removed"] = new_int_list_noise_removal


def get_mono_mass_peptide(df, noise=0.33, fixed_noise=5000):
    mz_cols = list(df.columns)
    np_arr = df.to_numpy()
    peptide_neutral_mass_list = []

    for row in np_arr:
        
        plain_peptide = row[mz_cols.index("plain_peptide")]
        modifications = row[mz_cols.index("modifications")]
        charge = int(row[mz_cols.index("charge_x")])
        
        total_mod_mass = get_mass_from_mods(modifications)
        peptide_neutral_mass = mass.calculate_mass(plain_peptide)+total_mod_mass

        peptide_neutral_mass_list.append(peptide_neutral_mass)

    df["peptide_neutral_mass"] = peptide_neutral_mass_list



def precursor_neutral_loss(df, precursor_charge_select_on, neutral_losses_dict,tol=10, intensity_cut=0.33, fixed_noise=5000): # we add nl_list that consists of neutral loss
    
    proton = mass.calculate_mass(formula='H+') #proron mono mass

    get_mono_mass_peptide(df, intensity_cut, fixed_noise)
    # we get intensity controled m/z and intensity and peptide neutral mass
    mz_cols = list(df.columns)
    np_arr = df.to_numpy()
    for nl in neutral_losses_dict.keys():
        nl_list = []
        for row in np_arr:
            peptide_neutral_mass = row[mz_cols.index("peptide_neutral_mass")]
            int_list = row[mz_cols.index("int_list_noise_removed")]
            mz_list = row[mz_cols.index("mz_list_noise_removed")]
            charge = int(row[mz_cols.index("charge_x")])

            if precursor_charge_select_on == "0":
                charge_init = charge
            else: 
                charge_init = 1

            nl_charge_list = []
            for x in range(charge_init, charge+1):

                # nl_mass = mass.calculate_mass(formula=nl.upper()) #computes nl mono mass
                nl_mass = neutral_losses_dict[nl]
                prec_mz_nl = (peptide_neutral_mass-nl_mass +(x*proton))/x #mz prec nl
                nl_charge_list.append(prec_mz_nl)

            match_nl = find_prec_match(nl_charge_list, mz_list, tol)
            nl_list.append(match_nl)
        df[nl+"_prec_loss"] = nl_list





# this section is used for server ...one ms2 per job
def ms2_DF(ms2File):
    x1 = ms2.IndexedMS2(ms2File,iterative=True)  #reading ms2 file using pyteomics iterative (bool, optional) – Defines whether iterative parsing should be used. It helps reduce memory usage at almost the same parsing speed. Default is True.
    ms2_df = pd.DataFrame([x for x in x1])  #dataframe of the mzXML file
    ms2_df["exp"] = os.path.basename(ms2File).split(".ms2")[0]
    ms2_df[["scan", "prec_MZ", "[M+H]+", "charge"]] = ms2_df.apply(parse_params_ms2, axis=1)
    ms2_df["spectrum"] = ms2_df["exp"]+"."+ms2_df["scan"].astype("str")+"."+ms2_df["charge"].astype("str")
    return ms2_df


# Precursor neutral loss functions 

def get_mass_from_mods(modifications):
    total_mod_mass = []
    each_residue_list = modifications.split(",")
    for x in each_residue_list:
        mod_val = float(x.split("_")[2])
        total_mod_mass.append(mod_val)
    return np.sum(total_mod_mass)

def noise_removal(mzlist, int_list, noise, fixed_noise): #33% as of https://pubs.acs.org/doi/full/10.1021/acs.jproteome.6b00487
    final_mz_list = []
    final_int_list = []
    
    baseline_intensity = np.max(int_list)
    noise = baseline_intensity*noise
    for i, mz in enumerate(mzlist):
        if (int_list[i] >= noise) & (int_list[i] >= fixed_noise):
            final_int_list.append(int_list[i])
            final_mz_list.append(mz)
    return final_int_list, final_mz_list
        

def find_prec_match(query, reference, tol):
    match = 0
    for mzlib in query:
        for index, masses in enumerate(reference):
            massshift = ppmCalc(float(mzlib), float(masses))
            if abs(massshift) < tol:
                match = 1
                break
    return match


def parse_params_ms2(row):
    params_dict = row.params
    scan = params_dict["scan"][0]
    prec_MZ = params_dict["precursor m/z"]
    neutral_mass = params_dict["neutral mass"][0]
    charge = int(params_dict["charge"][0])
    return pd.Series([scan, prec_MZ, neutral_mass, charge])


def ppmCalc(a, b, tol=10):
    a = float(a)
    b = float(b)
    #   massError = abs(a-b)*1e6/a
    massError = (b-a)*1e6/a  #calculates ppm error of the a(theoretical) from b(observed)
    return float(massError)


def neutral_loss_or_no_loss_prec(df):
    mz_cols = list(df.columns)
    np_arr = df.to_numpy()

    nl_cols = []
    for val in mz_cols:
        if "_prec_loss" in val:
            nl_cols.append(val)

    outcome_list_all = []

    for row in np_arr:
        outcome_list_rowwise = []
        for x in range(0, len(nl_cols)):
            outcome = row[mz_cols.index(nl_cols[x])]
            if outcome == "1":
                outcome_list_rowwise.append(nl_cols[x])

        outcome_list_all.append(",".join(outcome_list_rowwise))

    df["precursor_NL"] = outcome_list_all

def extract_cols_df(df, pattern_str="_prec_loss"):
    mz_cols = list(df.columns)
    nl_cols = []
    for cols in mz_cols:
        if pattern_str in cols:
            nl_cols.append(cols)
    return nl_cols

def precursor_nl_count(df,pattern_str="_prec_loss"):
    
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
            outcome_val = "No prec NL"
        else:
            outcome_val = "/".join(nl_row_list)
        
        store_nl_all.append(outcome_val)
    df["outcome_neutral_loss"] = store_nl_all


########### PLOTS
def make_pie_plot(df, tol, plot_axis="outcome",title = "Serine"):
    fig, ax = plt.subplots(figsize=(4,4))
    plt.style.use("ggplot")


    plt.rcParams['axes.edgecolor'] = "#010101"
    plt.rcParams['axes.facecolor'] = '#FFFFFF'


    plt.rcParams.update({'font.size': 10,'figure.max_open_warning': 0})
    ax = plt.axes(facecolor="w")

    y = np.array(df[plot_axis].value_counts().sort_index())
    mylabels = np.array(df[plot_axis].value_counts().sort_index().keys())

#     plt.title("Phosphopeptide containing modified {}".format(title))
    plt.title("{}".format(title))
    plt.pie(y, autopct='%1.0f%%', textprops={'color':"w",'weight':'bold', 'fontsize':10})

    plt.legend(labels = mylabels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    figurename = title+"_{}_precursor_neutral_loss.png".format(tol)
    
 
    plt.savefig(figurename, bbox_inches="tight", dpi=600 )