import sys, os
import pandas as pd
import glob
from os.path import dirname
import re
import numpy as np

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


def replace_symbol_with_val(df, stage_pep_dict, dynamic_mods_aa_dict):
    simplied_peptide_list = []
    df2 = df[["Peptides_original"]]
    mz_cols = list(df2.columns)
    np_arr = df2.to_numpy()
    for row in np_arr:
        original_pep = str(row[mz_cols.index("Peptides_original")])
        new_peptide = original_pep
        for stage_key in stage_pep_dict.keys():
            jump_modAA_dict, jump_mod_dict, sta_AA = getDynStatModsInfoPepXml(stage_pep_dict[stage_key])
            for delta_mass_key in  dynamic_mods_aa_dict.keys():
                add_mass = str(round(delta_mass_key, 2))
                for aa_residue in dynamic_mods_aa_dict[delta_mass_key]:
                    for symbols in jump_mod_dict.keys():
                        new_peptide = new_peptide.replace("{}{}".format(aa_residue, symbols) ,"{}[{}]".format(aa_residue, add_mass))
        simplied_peptide_list.append(new_peptide)
    df["Peptides_with_ptm_mass"] = simplied_peptide_list


def extract_spectrum_eval_tag(pep_xml, stage, tags_input_path):
    spectrum_list = []
    evalue_list = []
    tag_no_list = []
    stage_list = []
    mod_peptide = []
    spectrum = ""


    with open(pep_xml, "r") as f:
        for line in f:
            if '<spectrum_query spectrum=' in line.strip(): #spectrum found
                pattern = '<spectrum_query spectrum="(\w+\.\d+\.\d+\.\d+)"'
                spectrum=re.match(pattern, line.strip()).group(1)
                
               

                evalue = ""
                tag_num = ""
                mod_pep = ""
            if '<search_hit hit_rank="1"' in line:
                
                while "</search_hit>" not in line:  
                    #extract modified peptide info
                    if '<modification_info modified_peptide' in line.strip():
                        pattern_mod_pep = '<modification_info modified_peptide="(.+)" '
                        if mod_pep == "":
                            mod_pep=re.match(pattern_mod_pep, line.strip()).group(1)
                            mod_peptide.append(mod_pep)
                            spectrum_list.append(spectrum)
                            
                    if '<search_score name="expect"' in line.strip(): #evalue found
                        pattern_eval = '<search_score name="expect" value="(.+)"/>'
                        if evalue == "":
                            evalue=float(re.match(pattern_eval, line.strip()).group(1))
                            evalue_list.append(evalue)
    #                     print (evalue)

                    if '<search_score name="num_matched_tags"' in line.strip(): #evalue found
                        pattern_tag = '<search_score name="num_matched_tags" value="(.+)"/>'
                        if tag_num == "":
                            tag_num=float(re.match(pattern_tag, line.strip()).group(1))
                            tag_no_list.append(tag_num)

                    line = f.readline()
    stage_list = len(spectrum_list)*[stage]
    # print (len(spectrum_list), len(mod_peptide), len(evalue_list), len(tag_no_list), len(stage_list))

    if tags_input_path != "0":
        df = pd.DataFrame({"spectrum":spectrum_list,"mod_pep":mod_peptide,"expect":evalue_list, "no of matched tags":tag_no_list,"stage":stage_list})
    else:
        df = pd.DataFrame({"spectrum":spectrum_list,"mod_pep":mod_peptide,"expect":evalue_list, "stage":stage_list})
    
    df["spectrum__mod_pep"] = df.spectrum+"_"+df.mod_pep
    return df


def make_spectrum_single_jump_f(df):
    df["spectrum"] = df["Run#"]+"."+df["Scan#"].astype("str")+"."+df["z"].astype("str")


def make_pseudo_spectrum(row):
    spectrum = row.spectrum
    spec_split = spectrum.split(".")
    peudospec = spec_split[0]+"."+spec_split[1]+"."+spec_split[3]
    return peudospec


def quick_row_iterate(expDF):
    mz_cols = list(expDF.columns)
    np_arr = expDF.to_numpy()
    pseudo_spec_list = []
    for row in np_arr:
        spectrum = str(row[mz_cols.index("spectrum")])
        spec_split = spectrum.split(".")
        peudospec = spec_split[0]+"."+str(int(spec_split[1]))+"."+spec_split[3]
        pseudo_spec_list.append(peudospec)
    expDF["pseudo_spectrum"] = pseudo_spec_list



def merge_tag_files(tag_match_list): 
    #w010.10187.1.3_KIEESETIEDSSNQ[129]AAAR 
    df = pd.read_csv(tag_match_list[0], delimiter = "\t")
    stage=os.path.basename(dirname(dirname(dirname(tag_match_list[0]))))
    df["stage"] = stage
    df[["spectrum_no_pep","mod_pep"]] = df["spectrum"].str.split("_",expand=True)
    df[["Run#","Scan#","Precursor_rank","z"]] = df["spectrum_no_pep"].str.split(".",expand=True)
    df["spectrum_stage_key"]=df["Run#"]+"."+df["Scan#"].astype("str")+"."+df["z"].astype("str")+"__"+df["stage"]
    
   
    
    for x in range(1,len(tag_match_list)):
        df_2 = pd.read_csv(tag_match_list[x], delimiter = "\t")
        stage_2=os.path.basename(dirname(dirname(dirname(tag_match_list[x]))))
        df_2["stage"] = stage_2
        df_2[["spectrum_no_pep","mod_pep"]] = df_2["spectrum"].str.split("_",expand=True)
        df_2[["Run#","Scan#","Precursor_rank","z"]] = df_2["spectrum_no_pep"].str.split(".",expand=True)
        df_2["spectrum_stage_key"]=df_2["Run#"]+"."+df_2["Scan#"].astype("str")+"."+df_2["z"].astype("str")+"__"+df_2["stage"]
    
        frames = [df, df_2]
        df = pd.concat(frames)
    df.rename(columns = {"spectrum":"spectrum__mod_pep"}, inplace=True) #this is to merge in the evalue tag extracted file with same name key
    return df



def select_best_tag(df):
    mz_cols = list(df.columns)
    np_arr = df.to_numpy()
    best_tag_list = []
    best_score_list = []
    for row in np_arr:
        tags_all = row[mz_cols.index("Total_Tags")]
        tags_all = tags_all.replace("[","")
        tags_all = tags_all.replace("]","") #replace [and]
        tags_all = tags_all.replace(" ","") #replace whitespace
        tags_all = tags_all.replace("'","") #replace quotes
        tags_scores = row[mz_cols.index("Tag_scores_list")]
        tags_scores = tags_scores.replace("[","")
        tags_scores = tags_scores.replace("]","")
        
        tags_all_list = tags_all.split(",")
        tags_scores_list = tags_scores.split(",")
        tags_scores_list = [float(i) for i in tags_scores_list]
#         print (tags_scores_list)
        ind = np.argsort([-i for i in tags_scores_list])
        max_index = ind[0]
#         print (max_index)
        best_tag = tags_all_list[max_index]
        best_tag_score = tags_scores_list[max_index]
        
        best_tag_list.append(best_tag)
        best_score_list.append(best_tag_score)
    df["Top_scored_tag"] = best_tag_list
    df["Best_tag_scores"] = best_score_list
  