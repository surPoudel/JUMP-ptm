import pandas
import numpy as np
from itertools import product
import sys
import collections
import re

def reorder_psms(df): #this functions reorders psms based on weighted Evalues
    reordered_hits = []
    
#     start = time.time()
    mz_cols = list(df.columns)
    np_arr = df.to_numpy()
    for row in np_arr:
        rank_weightedEval_dict = {}
        # spectrum = str(row[mz_cols.index("spectrum")])
        spectrum = str(row[mz_cols.index("spectrum_jump")])
#         print (spectrum)
        search_hit = row[mz_cols.index("search_hit")]
        for index,val in enumerate(search_hit):
            weighted_e_value = val["search_score"]["WeightedEvalue"]
            rank_weightedEval_dict[index+1]= weighted_e_value
            
        list_of_scores = rank_weightedEval_dict.values()
        ind = np.argsort([-i for i in list_of_scores])
        search_hit_ranked_dict = {}
        
        for new_index, val_in in enumerate(ind):
            #replace with new rank
            search_hit[val_in]['hit_rank'] = new_index+1 
            search_hit_ranked_dict[new_index+1] = search_hit[val_in]
        
        od = collections.OrderedDict(sorted(search_hit_ranked_dict.items()))
        reordered_hits.append(list(od.values()))
    df["search_hit_reordered"] = reordered_hits
#     end = time.time()
#     print ("total time required to reorder search hits = ", end-start,"seconds\n")
     



#this function extracts dynamic and static modifications using pepxml files and stores them as the dictionary
def getPepXmlHeader(pepxml):
    header = []
    f = open(pepxml,"r") #read the file
    line = f.readline()
    while "</search_summary>" not in line: #end reading the file if this is sen
    # look for aminocaid modification keyword so that the modification infomration can be parsed
        header.append(line.rstrip())
        line = f.readline()
    header.append("</search_summary>")
    
    return header




def listToFile(pep_list, out):
    for x in pep_list:
        out.write(x+"\n")


#fucntion writePepXml changed to rank 1 and 2 
def writePepXml(pepxml,out,spectrum,prec_neutral_mass,assumed_charge,start_scan,end_scan,index,search_hit, spectrum_tag_total_dict):
    out.write('<spectrum_query spectrum="{}" start_scan="{}" end_scan="{}" precursor_neutral_mass="{}" assumed_charge="{}" index="{}">\n'.format(spectrum,start_scan,end_scan,prec_neutral_mass,assumed_charge,index) )
    out.write("  <search_result>\n")
    
    for hit in search_hit:
        rank = hit["hit_rank"]
        if (rank == 1) or (rank == 2):
            modified_peptide = hit["modified_peptide"]
            key_tag = spectrum+"_"+modified_peptide
            number_of_matchd_tags = 0
            try:
                number_of_matchd_tags = spectrum_tag_total_dict[key_tag]
            except:
                out.write("    No tag matched for {}".format(key_tag))
            out.write("   <search_hit ")
            out.write('{}="{}" '.format("hit_rank",hit["hit_rank"]))
            out.write('{}="{}" '.format("peptide",hit["peptide"]))
            hit_prot_dict = hit["proteins"][0]

            alternative_proteins = hit["proteins"][1:]

            #    <alternative_protein protein="sp|P52790|HXK3_HUMAN"/>


            out.write('{}="{}" '.format("peptide_prev_aa",hit_prot_dict["peptide_prev_aa"]))
            out.write('{}="{}" '.format("peptide_next_aa",hit_prot_dict["peptide_next_aa"]))
            out.write('{}="{}" '.format("protein",hit_prot_dict["protein"]))
            out.write('{}="{}" '.format("num_tot_proteins",hit["num_tot_proteins"]))
            out.write('{}="{}" '.format("num_matched_ions",hit["num_matched_ions"]))
            out.write('{}="{}" '.format("tot_num_ions",hit["tot_num_ions"]))
            out.write('{}="{}" '.format("calc_neutral_pep_mass",hit["calc_neutral_pep_mass"]))
            out.write('{}="{}" '.format("massdiff",hit["massdiff"]))
            out.write('{}="{}" '.format("num_tol_term",hit_prot_dict["num_tol_term"]))
            out.write('{}="{}" '.format("num_missed_cleavages",hit["num_missed_cleavages"]))
            out.write('{}="{}">\n'.format("num_matched_peptides",hit["num_matched_peptides"]))

            for alt_prot_dict in alternative_proteins:
                alt_protein = alt_prot_dict["protein"]
                out.write('    <alternative_protein protein="{}"/>\n'.format(alt_protein))

            out.write("    <modification_info ")
            out.write('{}="{}" '.format("modified_peptide",hit["modified_peptide"]))
            
            
            mod_amino_acids_list = hit["modifications"]
            for mod_amino_acids_dict in mod_amino_acids_list:
                
                if mod_amino_acids_dict["position"] != 0:
                    out.write('     <mod_aminoacid_mass position="{}" mass="{}"'.format(mod_amino_acids_dict["position"], mod_amino_acids_dict["mass"]))
                    if "static" in mod_amino_acids_dict.keys():
                        out.write(' static="{}"/>\n'.format(mod_amino_acids_dict["static"]))
                    if "variable" in mod_amino_acids_dict.keys():
                        out.write(' variable="{}" source="param"/>\n'.format(mod_amino_acids_dict["variable"]))
                    
                else:
                    out.write('{}="{}">\n'.format("mod_nterm_mass",mod_amino_acids_dict["mass"]))
            out.write('    </modification_info>\n')   
            search_score_dict = hit["search_score"]
            hit["search_score"]["num_matched_tags"] = number_of_matchd_tags

            for scoreKey in search_score_dict.keys():
                out.write('    <search_score name="{}" value="{}"/>\n'.format(scoreKey,search_score_dict[scoreKey]))
            
            out.write('   </search_hit>\n')

    

def makePepXML(df, pepxml,outputFile, spectrum_tag_total_dict, logFile):
    try:
        cmd = "rm {}".format(outputFile)
        os.system(cmd)
    except:
        # write_log (logFile,"File does not exist")
        print("File does not exist")
    
    header = getPepXmlHeader(pepxml)
    
    with open(outputFile, "a") as out:
        listToFile(header, out)
        
        mz_cols = list(df.columns)
        np_arr = df.to_numpy()
        for row in np_arr:
            # spectrum = str(row[mz_cols.index("spectrum")])
            spectrum = str(row[mz_cols.index("spectrum_jump")])
            prec_neutral_mass = float(row[mz_cols.index("precursor_neutral_mass")])
            assumed_charge = int(row[mz_cols.index("assumed_charge")])
            start_scan = int(row[mz_cols.index("start_scan")])
            end_scan = int(row[mz_cols.index("end_scan")])
            index = int(row[mz_cols.index("index")])

            search_hit =row[mz_cols.index("search_hit")]

            writePepXml(pepxml,out,spectrum,prec_neutral_mass,assumed_charge,start_scan,end_scan,index,search_hit,spectrum_tag_total_dict)
        out.write('  </search_result>\n')
        out.write(' </spectrum_query>\n')
        out.write(' </msms_run_summary>\n')
        out.write('</msms_pipeline_analysis>\n')  

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
                #print ("	",line.strip())
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



def keep_dynamic_mods_anyRank(dfMz, var_mass):
    dfMz = dfMz.dropna()
    
    spectrum_list = []
    mz_cols = list(dfMz.columns)
    np_arr = dfMz.to_numpy()
    for row in np_arr:
        spectrum = str(row[mz_cols.index("spectrum")])
        search_hits = row[mz_cols.index("search_hit")]

        for index,val in enumerate(search_hits):
            for x in val["modifications"]:
                if "variable" in x.keys():
                    for mod1 in var_mass:
                        if int(x["variable"]) == int(mod1):
                            spectrum_list.append(spectrum)
    spectrum_list = set(spectrum_list) #remove redundant spectrum
    dfMz_valid = dfMz.loc[dfMz.spectrum.isin(spectrum_list)]
    return dfMz_valid

#keep mods psms only

def keep_dynamic_mods(dfMz, var_mass):
    dfMz = dfMz.dropna()
    
    spectrum_list = []
    mz_cols = list(dfMz.columns)
    np_arr = dfMz.to_numpy()
    for row in np_arr:
        spectrum = str(row[mz_cols.index("spectrum")])
        search_hits = row[mz_cols.index("search_hit")]

        for index,val in enumerate(search_hits):
            if val["hit_rank"] == 1:
                for x in val["modifications"]:
                    if "variable" in x.keys():
                        for mod1 in var_mass:
                            if int(x["variable"]) == int(mod1):
                                spectrum_list.append(spectrum)
    spectrum_list = list(set(spectrum_list)) #remove redundant spectrum
    dfMz_valid = dfMz.loc[dfMz.spectrum.isin(spectrum_list)]
    return dfMz_valid,spectrum_list

def selectRanks(df):
#     start = time.time()
    all_hits = []
    mz_cols = list(df.columns)
    np_arr = df.to_numpy()
    for row in np_arr:
        
        search_hit = row[mz_cols.index("search_hit")]
        select_hits_list = []
        for index,val in enumerate(search_hit):
            if (val["hit_rank"] == 1) or (val["hit_rank"] == 2):
                select_hits_list.append(val)
        all_hits.append(select_hits_list)
    df["search_hit"] = all_hits
#     end = time.time()
#     print ("total time required to reorder search hits = ", end-start,"seconds\n")
     


def keep_dynamic_mods_anyRank(dfMz, var_mass):
    dfMz = dfMz.dropna()
    
    spectrum_list = []
    mz_cols = list(dfMz.columns)
    np_arr = dfMz.to_numpy()
    for row in np_arr:
        spectrum = str(row[mz_cols.index("spectrum")])
        search_hits = row[mz_cols.index("search_hit")]

        for index,val in enumerate(search_hits):
            for x in val["modifications"]:
                if "variable" in x.keys():
                    for mod1 in var_mass:
                        if int(x["variable"]) == int(mod1):
                            spectrum_list.append(spectrum)
    spectrum_list = set(spectrum_list) #remove redundant spectrum
    dfMz_valid = dfMz.loc[dfMz.spectrum.isin(spectrum_list)]
    return dfMz_valid

# def keep_dynamic_mods(dfMz, mod1, mod2):
#     dfMz = dfMz.dropna()
    
#     spectrum_list = []
#     mz_cols = list(dfMz.columns)
#     np_arr = dfMz.to_numpy()
#     for row in np_arr:
#         spectrum = str(row[mz_cols.index("spectrum")])
#         search_hits = row[mz_cols.index("search_hit")]

#         for index,val in enumerate(search_hits):
#             if val["hit_rank"] == 1:
#                 for x in val["modifications"]:
#                     if "variable" in x.keys():
#                         if (int(x["variable"]) == int(mod1)) or (int(x["variable"]) == int(mod2)):
#                             spectrum_list.append(spectrum)

#     dfMz_valid = dfMz.loc[dfMz.spectrum.isin(spectrum_list)]
#     return dfMz_valid


def write_log(logFile,*args):
    with open(logFile, "a") as log_file:
        line = ' '.join([str(a) for a in args])
        log_file.write(line+'\n')
        print(line)
