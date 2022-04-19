import os,sys
import numpy as np
from scipy import stats
from fragmenter import ionSeriesIonLossSpeRes,spectrumToDict
import pyteomics
from pyteomics import mass
from pepxml_generation_reorder import *

def make_mod_col(modifications_list):
    all_mods = []
    hydrogen = mass.calculate_mass(formula='H')
    for mod in modifications_list:
#         print (mod)
        if mod['position'] == 0:
            #n term static found
            all_mods.append("1_S_{}_n".format(mod["mass"]-hydrogen))
        elif ("variable" in mod.keys()) and ("static" in mod.keys()):
            all_mods.append("{}_V_{}_S_{}".format(mod["position"],mod["variable"], mod["static"]))
        
        
        elif "variable" in mod.keys():
            all_mods.append("{}_V_{}".format(mod["position"],mod["variable"]))
        elif "static" in mod.keys():
            all_mods.append("{}_S_{}".format(mod["position"],mod["static"]))
    return all_mods





#[m.start() for m in re.finditer('test', 'test test test test')]
#list_of_fragment_peptides has the fragment of peptides with the ptm in the middle
#all_matched_tags dictionary has key as spectrum, value = list of list of [tag, ion type]
def match_tag(df,scan_tag_dictionary, rank, spectrum, peptide, all_matched_tags,modified_peptide,min_tag_length,ionType="y+1"):
    
    #peptide = "".join(list(df.Seq)) #gets the sequence from the fragmenter
    for key in scan_tag_dictionary[spectrum].keys():
        #initiate tag_score, rtflank_tag
#         print (key)
        if (key in peptide) or (key[::-1] in peptide):
            #get tag file data for that key 
            tag_file_data = scan_tag_dictionary[spectrum][key]#gets a list of Evalue, rtflanking mass and rank of the tag
            # print (tag_file_data)
            tag_score = tag_file_data[0]
            rtflank_tag = tag_file_data[1]
            rtflank_tag_rank = tag_file_data[2]

#             index_val_list = [list(m.span()) for m in re.finditer(key, peptide)][0]
            index_val_list = [i for i in range(len(peptide)) if peptide.startswith(key, i)]
            index_val_list2 = [i for i in range(len(peptide)) if peptide.startswith(key[::-1], i)]
            
            index_val_list+=index_val_list2
                

#             print (key, index_val_list)
#             print (key+"\t",index_val_list)
#             index_val_list = [m.start() for m in re.finditer(key, peptide)]
#             index_val = peptide.index(key) #gets the index of start of tag
#             print (index_val_list)
            
            for index_val in index_val_list:
                if ionType == "b+1":
                    rightSideMass_theo = 0
                    if index_val+len(key) < df.shape[0]:

                        tag_index = index_val+len(key)-1
                        start_position = index_val+1
                        end_position = start_position+len(key)-1

                        rightSideMass_theo = df.iloc[tag_index][ionType]
#                     print (rightSideMass_theo)
                else:
                    tag_index = index_val
                    start_position = index_val+1
                    end_position = start_position+len(key)-1

                    rightSideMass_theo = df.iloc[tag_index][ionType]
#                     print (rightSideMass_theo)
                if rightSideMass_theo != '---':
                    if int(float(rtflank_tag)) == int(float(rightSideMass_theo)):
                        #match found #update tag dictionary
                        
                        if len(key) >=min_tag_length:
                            # print ("    We have a tag match key = {} value = {}".format(spectrum+"_"+modified_peptide,[key,ionType,position,tag_score,rtflank_tag_rank]))

                            valAddKey(all_matched_tags, spectrum+"_"+modified_peptide, [key,ionType,rank,tag_index, start_position, end_position,tag_score,rtflank_tag_rank])
                            




def dividePeptideFromMods(modifications, peptide, delta_mass = 79.966331):
    list_of_fragment_peptides = []
    peptide_fragments = [] #index 0 is before ptm (left part) and index 1 is with ptm (right part)
    for val in modifications:
        if 'variable' in val.keys():
            if val['variable'] == delta_mass: #phospho_found
                position = val["position"] #position of ptm
                #split the peptide to 2 parts 
#                 peptide_fragment = [peptide[0:position],peptide[position:]]
                peptide_length = len(peptide)
                add_len = peptide_length - position
                list_of_fragment_peptides.append(peptide[0:position]+("Z"*add_len))
                list_of_fragment_peptides.append("Z"*(peptide_length-add_len)+peptide[position:])
#                 list_of_fragment_peptides.append(peptide_fragment)
    try:
        peptide_fragments=[list_of_fragment_peptides[0], list_of_fragment_peptides[-1]]
    except:
        print ("No modifications in the peptide")
    return  peptide_fragments#two extreme peptide sequence fragments with ptms in middle




def makedir(dir1):
  cmd = "mkdir "+dir1
  try:
    os.system(cmd)
  except:
    print ("{} Directory exits!!".format(dir1))





def ppmCalc(a, b, tol=10):
    a = float(a)
    b = float(b)
    #   massError = abs(a-b)*1e6/a
    massError = (b-a)*1e6/a  #calculates ppm error of the a(theoretical) from b(observed)
    return float(massError)


def match_precursor(jump_comet_link_dict, comet_key):
    spectrum = ""
    for key in jump_comet_link_dict.keys():
        comet_key_split_under = comet_key.split("__")
        comet_key_split = comet_key_split_under[0].split(".")
        key_phrase = comet_key_split[0]+"."+comet_key_split[1]+"."+comet_key_split[2]+"."+comet_key_split[3]

        if key_phrase in key:
            # print (key_phrase)
            # print (jump_comet_link_dict[key])
            #potential match found
            #get comet key neutral mass
            key_split = key.split("_")
            neutral_mass_tag = float(key_split[-1])
            neutral_mass_comet = float(comet_key_split_under[-1])
            # print (neutral_mass_tag,"\t",neutral_mass_comet)
            massdiff = ppmCalc(neutral_mass_tag, neutral_mass_comet)
            # print (massdiff)
            if abs(massdiff) < 50:
                # print (massdiff, key, jump_comet_link_dict[key])
                spectrum = jump_comet_link_dict[key]
                break

    return spectrum




def get_tag_matched(df, scan_tag_dictionary,jump_comet_link,all_matched_tags, logFile, min_tag_length,dataType = "whole"): #df is one row dataframe 
#     start = time.time()
    #to store new outfiles missing in tags file
    arbitrary_outfile = []
    cnt = 0
    jump_spectrum_list = []

    mz_cols = list(df.columns)
    np_arr = df.to_numpy()
    for row in np_arr:
        comet_spectrum = str(row[mz_cols.index("spectrum")])

        if comet_spectrum in arbitrary_outfile:
            cnt=+1
        else:
            arbitrary_outfile.append(comet_spectrum)
            cnt=1
        neutral_mass = str(row[mz_cols.index("precursor_neutral_mass")])   
        charge = int(row[mz_cols.index("assumed_charge")]) 

        MH = float(neutral_mass)+1.00782503207 #mass of hydrogen
        comet_key = comet_spectrum+"__"+str(MH)
        search_hit = row[mz_cols.index("search_hit")]
        jump_key = match_precursor(jump_comet_link, comet_key)
        
        # if comet_spectrum == "w005.14658.14658.2":
        #     print ("comet key = " ,comet_key)
        #     print ("jump key = ",jump_key)

        #updates jump key if tags do not have the scan
        if jump_key == "":
            # print("    {} does not exist in tag file. So, arbitrary jump like outfile created".format(comet_key))
            write_log(logFile,"    {} does not exist in tag file. So, arbitrary jump like outfile created".format(comet_key))
            comet_split = comet_spectrum.split(".")
            jump_key = "{}.{}.{}.{}".format(comet_split[0],comet_split[1],cnt,charge)
            
            # if comet_spectrum == "w005.14658.14658.2":
            #     print ("comet key = " ,comet_key)
            #     print ("jump key = ",jump_key)
                

        # print ("The jump key for {} is {}".format(comet_spectrum, jump_key))
        jump_spectrum_list.append(jump_key) #replace comet spectrum with jump spectrum

        for index,val in enumerate(search_hit):
            rank = int(val["hit_rank"])
            if (rank == 1) or (rank == 2):
                modified_peptide = val["modified_peptide"]
                # comet_evalue = get_comet_scores("expect", val["search_score"])
    #             comet_evalue = get_comet_scores("xcorr", val["search_score"])
        

                if jump_key in scan_tag_dictionary.keys():
    #                 print (spectrum)
    #                 print (comet_evalue)
                    modifications = val["modifications"] #list of dictionary
                    all_mods = make_mod_col(modifications)
                    peptide = val["peptide"]
                    modifications_combine = ",".join(all_mods)
                    massPosDict1 = spectrumToDict(modifications_combine)
                    df_pep = ionSeriesIonLossSpeRes(peptide,maxcharge=2,massPosDict=massPosDict1,useMod ="Yes")
                    
                    
                    if dataType != "whole":
                        #this is required for the ptms to get the fragment from the mods
                        list_of_fragment_peptides = dividePeptideFromMods(modifications, peptide, 79.966331)
                    else:
                        list_of_fragment_peptides = [peptide, peptide] #just to make the compatible
    #                 print (list_of_fragment_peptides)
                    if len(list_of_fragment_peptides)!=0:
                        #left flanks scores and information
                        
                        match_tag(df_pep,scan_tag_dictionary, rank, jump_key, list_of_fragment_peptides[0],all_matched_tags,modified_peptide,min_tag_length,"b+1")
                        match_tag(df_pep,scan_tag_dictionary, rank, jump_key, list_of_fragment_peptides[-1],all_matched_tags,modified_peptide,min_tag_length)
                    
    return jump_spectrum_list
    