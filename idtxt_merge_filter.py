import sys, os
import pandas as pd
import glob

import argparse
import numpy as np
import re
from logFunctions import *
from stagewise_comet_search import *
import subprocess
import pyteomics as pyteomics
from pyteomics import mass

from os.path import dirname


#This fucntion computes modifications with all the static and dynamic represenatation
def computeModifications(row, jump_mod_dict, sta_AA, peptide = "Peptides"):
    mod_peptide = row[peptide] #Peptide sequence predicted by the software
    # print (mod_peptide,"\t",row.Outfile)
    # print (row.Outfile)
    peptideNoFlank = mod_peptide.split(".")[1] #Flanking peptide
    pepNoFlankList = list(peptideNoFlank) #flanking peptide is converted to the list, so the symbols are also member of list along with amino acids

    plain_peptide_list =[] #iniiates plain peptide formation
    for aa in pepNoFlankList:
        if aa in pyteomics.parser.std_amino_acids: #this looks at the standard amino acids in pyteomics and if the value in list is not amino acid for example * or other symbols, then it will discard those and just adds amino acdis
            plain_peptide_list.append(aa) 
    plain_peptide = "".join(plain_peptide_list) #creates a new string or amino acids
    # print (plain_peptide)
    dynAA_count = 0 #counts the number of dynamic modifcaiton in the peptide sequence
    for key in jump_mod_dict.keys(): 
        if key in mod_peptide:
            dynAA_count +=1 #updates dynamic modification if it finds the symbol in modified peptide
    
    # print (dynAA_count)
    
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
    # print (mods)
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




#this function is different than validator fucntion as the masses are as list here and not added
def spectrumToDict(spectrum):
    dict1 = {}
    spectrumCommaSplit = spectrum.split(",")
    for x in spectrumCommaSplit:
        y=x.split("_")
        if y[0] not in dict1.keys():
            dict1[y[0]] = [str(y[2])]
        else:
            dict1[y[0]].append(str(y[2]))
    dict2 = {}

    for key in dict1.keys():
        value = "+".join(list(dict1[key])) #for fully tryptic peptides
        # value = "+".join(list(set(dict1[key])))
        dict2[key] = value
    return dict2


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



def addModValInPepSeq(plain_peptide,massPosDict): #massPosDict example is massPosList[0]
    modified_peptide = []
    for index, aa in enumerate(plain_peptide):
        pos = str(index+1)
        if pos in massPosDict.keys():
            massShift = massPosDict[pos]
            aa2 = aa+"("+str(massShift)+")"
            modified_peptide.append(aa2)
        else:
            modified_peptide.append(aa)
    return ("".join(modified_peptide))


#jump -q only looks at the peptide sequence so if same peptide exists from different stages, it still considers same peptide and produces same quantification
def make_unique_peptide(df, pep_col = "Peptides"):
    all_peptides = list(df[pep_col])
    non_redun_peptides = list(set(all_peptides))
    replica_all_peptides = all_peptides
    
    mod_peptides = []
    for index,pep in enumerate(non_redun_peptides):

        indices = [i for i, x in enumerate(all_peptides) if x == pep]
        pep_split = pep.split(".")
        for mini_index, val in enumerate(indices):
            pep_combine = "{}.{}{}.{}".format(pep_split[0],pep_split[1],(mini_index)*"#",pep_split[2])
            replica_all_peptides.pop(val)
            replica_all_peptides.insert(val, pep_combine)
    
    df["Peptides_unique_made"] = replica_all_peptides



def make_unique_peptide_addZ(row, peptide="Peptides"):
    peptide = row[peptide]
    stage = row.ptm_stage
    add_z = int(stage.split("_")[-1])
    pep_split = peptide.split(".")
    pep_combine = "{}.{}{}.{}".format(pep_split[0],pep_split[1],add_z*"Z",pep_split[2])
    return pep_combine


def merge_publication_jump_f_tables(id_files, out_Folder, stage0_idtxt,keep_spectrum):
    #to merge the stage0 publication tables
    #Run#    Scan#   m/z     z
    if "id_uni_pep" in id_files[0]:
        stage_0_pep = dirname(stage0_idtxt)+"/publications/id_uni_pep.txt"
    else:
        stage_0_pep = dirname(stage0_idtxt)+"/publications/id_all_pep.txt"

    stage_0pepDF = pd.read_csv(stage_0_pep, delimiter="\t", skiprows=return_skiprows(stage_0_pep, "\t", "Peptides"))
    stage_0pepDF["spectrum"] = stage_0pepDF["Run#"]+"."+stage_0pepDF["Scan#"].astype("str")+"."+stage_0pepDF["z"].astype("str")
    stage_0pepDF["ptm_stage"] = os.path.basename(dirname(dirname(stage0_idtxt)))
    stage_0pepDF["key"] = stage_0pepDF.ptm_stage+"_"+stage_0pepDF.spectrum


    pepDF = pd.read_csv(id_files[0], delimiter="\t", skiprows=return_skiprows(id_files[0], "\t", "Peptides"))
    pepDF["spectrum"] = pepDF["Run#"]+"."+pepDF["Scan#"].astype("str")+"."+pepDF["z"].astype("str")
    pepDF["ptm_stage"] = os.path.basename(dirname(dirname(dirname(id_files[0]))))
    pepDF["key"] = pepDF.ptm_stage+"_"+pepDF.spectrum

    pepDF = pepDF.loc[pepDF["key"].isin(keep_spectrum)]

    all_cols = list(pepDF.columns)+["Peptides_original"]


    # #pepDF["Peptides"] = pepDF["Peptides"]+"_"+pepDF["ptm_stage"]
    pepDF["Peptides_unique_made"] = pepDF.apply(make_unique_peptide_addZ, peptide="Peptides", axis=1) 
    # make_unique_peptide(pepDF) #this creates a unique peptide name with # at the end depending on number of redundant peptide. First will have same mods, second will have 1 # and third will have 2# at the end and so on
    # #a new column called Peptides_unique_made forms
    pepDF.rename(columns={"Peptides":"Peptides_original"}, inplace=True)
    pepDF.rename(columns={"Peptides_unique_made":"Peptides"}, inplace=True)

    pepDF = pepDF[all_cols]

 

    for pep_file in id_files[1:]:
        pepDF_1 = pd.read_csv(pep_file, delimiter="\t", skiprows=return_skiprows(pep_file, "\t", "Peptides"))
        pepDF_1["spectrum"] = pepDF_1["Run#"]+"."+pepDF_1["Scan#"].astype("str")+"."+pepDF_1["z"].astype("str")
        pepDF_1["ptm_stage"] = os.path.basename(dirname(dirname(dirname(pep_file))))
        pepDF_1["key"] = pepDF_1.ptm_stage+"_"+pepDF_1.spectrum
        pepDF_1 = pepDF_1.loc[pepDF_1["key"].isin(keep_spectrum)]

        # pepDF_1["ptm_stage"] = os.path.basename(dirname(dirname(dirname(pep_file))))
        #pepDF_1["Peptides"] = pepDF_1["Peptides"]+"_"+pepDF_1["ptm_stage"]

        #make_unique_peptide(pepDF_1) #this creates a unique peptide name with # at the end depending on number of redundant peptide. First will have same mods, second will have 1 # and third will have 2# at the end and so on
        #a new column called Peptides_unique_made forms
        pepDF_1["Peptides_unique_made"] = pepDF_1.apply(make_unique_peptide_addZ, peptide="Peptides", axis=1) 

        pepDF_1.rename(columns={"Peptides":"Peptides_original"}, inplace=True)
        pepDF_1.rename(columns={"Peptides_unique_made":"Peptides"}, inplace=True)

        pepDF_1 = pepDF_1[all_cols]


        frames = [pepDF, pepDF_1]
        pepDF = pd.concat(frames)
    #add stage 0 too in final pepDF

    # pepDF2 = pepDF.loc[pepDF["key"].isin(keep_spectrum)]

    #replace Jscore with Xcorr  and dJn  with dCn for stage 0 id files
    stage_0pepDF = stage_0pepDF.rename(columns={"Jscore":"Xcorr","dJn":"dCn"})

    frames = [pepDF, stage_0pepDF]
    pepDF = pd.concat(frames)
    ####

    makedirectory(out_Folder+"/publications")
    outputFileName = os.path.basename(id_files[0])

    pepDF.to_csv(out_Folder+"/publications/"+outputFileName, sep="\t", index=None)
    
    
def combine_jump_f_stages_modPeptides(all_stage_idtxts, stage0_idtxt):
    print ("Merging {}".format(all_stage_idtxts[0]))

    f1 = open(all_stage_idtxts[0],"r")
    head1=f1.readline()

    pepxml = glob.glob(dirname(dirname(all_stage_idtxts[0]))+"/*/*.pep*")[0]
#     print (pepxml)
    jump_modAA_dict, jump_mod_dict, sta_AA = getDynStatModsInfoPepXml(pepxml)

    # print (jump_modAA_dict, jump_mod_dict, sta_AA)

    idtxtdf = pd.read_csv(all_stage_idtxts[0], delimiter=";", skiprows=return_skiprows(all_stage_idtxts[0], ";", "Peptide"))
    
    original_cols = list(idtxtdf.columns)
    
    #keep the ptms only .. remove met oxidation if given
    symbols = list(jump_mod_dict.keys())
#     print (symbols)
    if "*" in symbols:
        symbols.remove("*") #removes methionine oxidation symbol

    ptms_symbols = "|".join(symbols)
#     print (ptms_symbols)
    idtxtdf = idtxtdf.loc[idtxtdf["Peptide"].str.contains(ptms_symbols)]
#     print (idtxtdf.shape)
#     idtxtdf["spectrum"] = idtxtdf.apply(createOutfile, df=idtxtdf, axis=1)

#     idtxtdf[["exp","scan","charge"]] = idtxtdf["spectrum"].str.split(".",expand=True)

    idtxtdf["Peptides"] = idtxtdf["Peptide"]

    # print (idtxtdf.columns)

    #For modification or PTMs information on the consensus library
    idtxtdf[["plain_peptide","modifications"]] = idtxtdf.apply(computeModifications, jump_mod_dict=jump_mod_dict,sta_AA=sta_AA,axis=1)
    idtxtdf["massPosDict"] = idtxtdf.modifications.apply(lambda x: spectrumToDict(x))
    idtxtdf["modified_peptide_with_mass"] = idtxtdf.apply(lambda x: addModValInPepSeq(x.plain_peptide,x.massPosDict), axis=1)
    idtxtdf["ptm_stage"] = os.path.basename(dirname(dirname(all_stage_idtxts[0])))

    #idtxtdf["Peptide"] = idtxtdf["Peptide"]+"_"+idtxtdf["ptm_stage"]
    #make_unique_peptide(idtxtdf,"Peptide") #this creates a unique peptide name with # at the end depending on number of redundant peptide. First will have same mods, second will have 1 # and third will have 2# at the end and so on
    #a new column called Peptides_unique_made forms
    idtxtdf["Peptides_unique_made"] = idtxtdf.apply(make_unique_peptide_addZ, peptide="Peptide", axis=1) 


    idtxtdf.rename(columns={"Peptide":"Peptides_original"}, inplace=True)
    idtxtdf.rename(columns={"Peptides_unique_made":"Peptide"}, inplace=True)


    for idtxtFile in all_stage_idtxts[1:]:
        print ("Merging {}".format(idtxtFile))
        pepxml = glob.glob(dirname(dirname(idtxtFile))+"/*/*.pep*")[0]
        jump_modAA_dict_1, jump_mod_dict_1, sta_AA_1 = getDynStatModsInfoPepXml(pepxml)
        # print (jump_modAA_dict_1, jump_mod_dict_1, sta_AA_1)
        idtxtdf_1 = pd.read_csv(idtxtFile, delimiter=";", skiprows=return_skiprows(idtxtFile, ";", "Peptide"))

        #keep the ptms only .. remove met oxidation if given
        symbols_1 = list(jump_mod_dict_1.keys())
        if "*" in symbols_1:
            symbols_1.remove("*") #removes methionine oxidation symbol
    
        ptms_symbols_1 = "|".join(symbols_1)
        # print (ptms_symbols_1)
        idtxtdf_1 = idtxtdf_1.loc[idtxtdf_1["Peptide"].str.contains(ptms_symbols_1)]

#         idtxtdf_1["spectrum"] = idtxtdf_1.apply(createOutfile, df=idtxtdf_1, axis=1)

#         idtxtdf_1[["exp","scan","charge"]] = idtxtdf_1["spectrum"].str.split(".",expand=True)
        idtxtdf_1["Peptides"] = idtxtdf_1["Peptide"]
        # print (idtxtdf_1.shape)
        # print (idtxtdf_1.columns)
        #For modification or PTMs information on the consensus library
        idtxtdf_1[["plain_peptide","modifications"]] = idtxtdf_1.apply(computeModifications, jump_mod_dict=jump_mod_dict_1,sta_AA=sta_AA_1,axis=1)
        idtxtdf_1["massPosDict"] = idtxtdf_1.modifications.apply(lambda x: spectrumToDict(x))
        idtxtdf_1["modified_peptide_with_mass"] = idtxtdf_1.apply(lambda x: addModValInPepSeq(x.plain_peptide,x.massPosDict), axis=1)
        idtxtdf_1["ptm_stage"] = os.path.basename(dirname(dirname(idtxtFile)))
        #idtxtdf_1["Peptide"] = idtxtdf_1["Peptide"]+"_"+idtxtdf_1["ptm_stage"]

        idtxtdf_1["Peptides_unique_made"] = idtxtdf_1.apply(make_unique_peptide_addZ, peptide="Peptide", axis=1)  #this creates a unique peptide name with # at the end depending on number of redundant peptide. First will have same mods, second will have 1 # and third will have 2# at the end and so on
        #a new column called Peptides_unique_made forms
        idtxtdf_1.rename(columns={"Peptide":"Peptides_original"}, inplace=True)
        idtxtdf_1.rename(columns={"Peptides_unique_made":"Peptide"}, inplace=True)


        frames = [idtxtdf, idtxtdf_1]
        idtxtdf = pd.concat(frames)
        
    final_cols = original_cols+["plain_peptide","modifications","modified_peptide_with_mass","ptm_stage","Peptides_original"]
    #stag0 idtxt merge for loading bias correction for jump -q
    
    print ("Merging {}".format(stage0_idtxt))
    idDF = pd.read_csv(stage0_idtxt, delimiter=";", skiprows=return_skiprows(stage0_idtxt, ";", "Peptide"))
    idDF["ptm_stage"] = "Stage_0"
    idtxtdf = idtxtdf[final_cols].append(idDF)

    return head1,idtxtdf[final_cols]

#program to select the strategy to combine filter result
#Options for QC after concatenation of jump -f stagewise results: 1: Xcorr based priority (best xcorr scan retains) [DEFAULT], 2: Stage priority (Stage 1 wins over Stage2), 3: Keep all result (no removal)- Chances few scans could have different peptides but identical quantification), 4: Remove all redundant scans
#combine_filter_result = 1


def get_non_redundant_spectra_xcorr(df2):

    df = df2.copy()

    df["spectrum"] = df.apply(createOutfile, df=df, axis=1)
    df["XCorr"] = df.XCorr.astype("float")
    df = df.sort_values(by=["XCorr"], ascending=False)
    df=df.drop_duplicates(subset=["spectrum"], keep="first")

    df["key"] = df.ptm_stage+"_"+df.spectrum

    keep_spectrum = list(set(list(df["key"])))

    # print (keep_spectrum)
    return keep_spectrum


def get_non_redundant_spectra_stagewise(df):
    df["spectrum"] = df.Outfile.apply(lambda x: os.path.basename(x))
    
    df = df.sort_values(by=["spectrum","ptm_stage"], ascending=[True,True])
    


    df["spectrum"] = df.apply(createOutfile, df=df, axis=1)
    df2=df.drop_duplicates(subset=["spectrum"], keep="first")
    df2["key"] = df.ptm_stage+"_"+df.spectrum
    
    #unique keys 
    uniqueKeys = list(df2.key.drop_duplicates())
    #subset the df for the output
    df = df.loc[df.key.isin(uniqueKeys)]
    keep_spectrum = list(df["key"])
    return keep_spectrum




def get_non_redundant_spectra_remove_all_duplicates(df):

    df["spectrum"] = df.Outfile.apply(lambda x: os.path.basename(x))
    df["key"] = df.ptm_stage+"_"+df.spectrum
 
    
    #check df2 to remove duplicate spectrum
    uniqueSpectrum = list(df.spectrum.drop_duplicates())
    
    #subset the df for the output
    df = df.loc[df.spectrum.isin(uniqueSpectrum)]
    
    keep_spectrum = list(df["key"])
    return keep_spectrum


def combine_filter_result_consolidate(df, keep_spectrum):#input is the dataframe df that is the output from combine_jump_f_stages_modPeptides(all_stage_idtxts)
    reqd_cols = df.columns
    # print(reqd_cols)
    
    #df["spectrum"] = df.Outfile.apply(lambda x: os.path.basename(x))
    
    #w072.18564.1.3.out
    df["spectrum"] = df.apply(createOutfile, df=df, axis=1)

    df["key"] = df.ptm_stage+"_"+df.spectrum
    # df2=df.drop_duplicates(subset=["spectrum"], keep="first")
    
    df = df.loc[df["key"].isin(keep_spectrum)]


    return df[reqd_cols]

