import pyteomics
from pyteomics import pepxml
import pandas as pd
import numpy as np
import os,sys

from fragmenter import *
from pepxml_generation_reorder import *
from tagExtraction import *
from target_decoy_qc_plots import *
from assign_tag_combine_score import *
import configparser
config = configparser.ConfigParser()

import time
import pickle

import os.path
from os.path import dirname


# pep_xml = "/home/spoudel1/PanPTM_Paper_2021/Tag_based_QC/Same_DB_comet_search/w010/w010.pep.xml"
# tag_file = "/home/spoudel1/PanPTM_Paper_2021/w010/w010.1.tags"
# test_tag = "/home/spoudel1/PanPTM_Paper_2021/w010/w010_9741.tags"

inpath = os.getcwd()
params_file = sys.argv[1]

config.read(params_file)

pep_xml = config["tagqc"]["pep_xml"]
tag_file = config["tagqc"]["tag_file"]
start_scan = config["tagqc"]["start_scan"]
end_scan = config["tagqc"]["end_scan"]
min_tag_length = config["tagqc"]["min_tag_length"]

logFile = inpath+"/tag_match.log"

# pep_xml = sys.argv[1] #pep_xml file example: /home/spoudel1/PanPTM_Paper_2021/Tag_based_QC/Same_DB_comet_search/w010/w010.pep.xml
# tag_file = sys.argv[2] #tag_file: /home/spoudel1/PanPTM_Paper_2021/w010/w010.1.tags
# start_scan = 20000
# end_scan = 30000

'''
#extract positive control set of psms (spectrum)
positive = "Top_100_psms_jump_search.xlsx"
positiveDF = pd.read_excel(positive)
positveSpectrum = list(positiveDF.spectrum)
'''

jump_modAA_dict, jump_mod_dict, sta_AA = getDynStatModsInfoPepXml(pep_xml)
# mod1= 79.966331
# mod2 = 81.5000

# print (jump_modAA_dict)
# print (jump_mod_dict)
# print (sta_AA)

var_mass = list(jump_modAA_dict.keys())
# var_mass = [79.966331]
# print (var_mass)

fraction_name = os.path.basename(pep_xml).split(".")[0]

os.chdir(os.getcwd())
out_dir = "Results_start_scan_{}_end_scan_{}_min_tag_len_{}".format(start_scan,end_scan,min_tag_length)
makedir(out_dir)

output_directory = dirname(pep_xml)

write_log (logFile,"    Reading pepxml file and storing as dataframe")
# print ("    Reading pepxml file and storing as dataframe")
x1 = pyteomics.pepxml.read(pep_xml)  #reading pepxml file using pyteomics
dfMz = pd.DataFrame([x for x in x1])  #dataframe of the pepxml file
write_log (logFile,".....Reading Done\n")

# dfMz_valid_test2 = dfMz.loc[dfMz.spectrum.isin(["w005.14658.14658.2"])]
# print (dfMz_valid_test2)

#remove isna search hits

#here keep only ptms if the data search is ptm type
# dfMz_valid = keep_dynamic_mods(dfMz, mod1, mod2)
# dfMz_valid = keep_dynamic_mods(dfMz, var_mass)

write_log (logFile,"Keeping Rank1 psms with ptms only")
write_log (logFile,"    Initial spectrum # {}".format(dfMz.shape[0]))
dfMz_valid,spectrum_list = keep_dynamic_mods(dfMz, var_mass)
# print (spectrum_list)

#selectRanks(dfMz_valid)
# dfMz_valid_test = dfMz_valid.loc[dfMz_valid.spectrum=="w005.14658.14658.2"]
# print (dfMz_valid_test)
write_log (logFile,"Storing tag information and generating outfile links from comet to jump")
#scan tag has jump spectrum linked to tags, jump_comet_link has comet spectrum + MH link to jump tag
# scan_tag, tag_score_list, jump_comet_link = scan_tag_dictionary(tag_file)

scan_tag, tag_score_list, jump_comet_link = scan_tag_dictionary(tag_file)


write_log (logFile,"    Tag information stored")

#select scan_tag and jump_comet_link
#scan_tag_dict, tag_score_list, jump_comet_link_dict

# for key in spectrum_list:
#     jump_outfile = jump_comet_link_dict[key]
#     jump_comet_link[key] = jump_outfile
#     scan_tag[jump_outfile]  = scan_tag_dict[jump_outfile]


tag_mean = np.mean(tag_score_list)
tag_std = np.std(tag_score_list)


write_log (logFile,"    THe mean of tag score of all tags is {}".format(tag_mean))
write_log (logFile,"    THe standard deviation of tag score of all tags is {}".format(tag_std))



write_log (logFile,"    ptm containing rank1 psms retained")
write_log (logFile,"    Total retained spectrum # {}".format(dfMz_valid.shape[0]))

# print (dfMz_valid)

#get scans to test


# print ("    Check mean result {}\n\tstd dev {}\n\n\n".format(comet_eval_mean,comet_eval_std))


if end_scan != "0": 
    dfMz_valid = dfMz_valid.iloc[int(start_scan):int(end_scan)]

# print (dfMz_valid)


'''
#top 100 test set
dfMz_valid = dfMz_valid.loc[dfMz_valid.spectrum.isin(positveSpectrum)]
single spectrum
dfMz_valid=dfMz_valid.loc[dfMz_valid.spectrum == "w010.09935.09935.2"]
'''

os.chdir(out_dir)

# file_Check ="{}_reordered_final.pickle".format(fraction_name) 
# if os.path.exists(file_Check):
#     print ("The program found already evaluated {}. This file will be used to make dataframe with combined score".format(file_Check))
#     with open(file_Check, 'rb') as handle:
#         dfMz_valid = pickle.load(handle)

# else:
start = time.time()
df_split = np.array_split(dfMz_valid, 10)
cnt = 1
n=1
all_matched_tags = {}

jump_spectrum_all = []

write_log (logFile,"\nExecuting Tag Match Program")
for df_x in df_split:
    #function get_tag_matched  is changed to rank 1 and rank2 
    jump_spectrum_list = get_tag_matched(df_x, scan_tag,jump_comet_link,all_matched_tags,logFile,int(min_tag_length))
    write_log (logFile,"    total spectrums matched {} ".format(len(jump_spectrum_list)))
    jump_spectrum_all += jump_spectrum_list
    cnt+=1
    
    if cnt == n*20:
        write_log (logFile,"    total frames {} ".format(cnt))
        n+=1
end = time.time()
write_log (logFile,"total time required to reorder search hits = ", end-start,"seconds\n")

dfMz_valid["spectrum_jump"] = jump_spectrum_all


with open('{}_reordered_final.pickle'.format(fraction_name), 'wb') as handle:
    pickle.dump(dfMz_valid, handle, protocol=pickle.HIGHEST_PROTOCOL)

#get all matched tags
#all_matched_tags, spectrum, [key+"\t"+ionType,"\t",position]
# print (all_matched_tags)
write_log (logFile,"    Reporting tags per spectrum in the tag output file")
out_file = "spectrum_tag_count.txt"


report_tag(all_matched_tags, out_file)

#use out_file to compute total matched tags number
#If a tag starts from same position and has a better score and another tag is shorter also starts from same position
#the tag with longer position is retained .. thus total tags number are computed. A subset tag is emoved

#generate dataframe using out_file 

def my_agg(x):
    names = {
        'Total_Tags': x['tag'].to_list(),
        'Total_Tag#': x['tag'].count(),
        'Tag_ions_list': x['ionType'].to_list(),
        'Tag_starts_list': x['tag_start_position'].to_list(),
        'Tag_ends_list': x['tag_end_position'].to_list(),
        'Tag_scores_list': x['tag_score'].to_list(),
        'Tag_ranks_list': x['tag_rank'].to_list(),
        }
    return pd.Series(names)

def total_tags_compute(spectagfile):
    df_tag_spec = pd.read_csv(spectagfile, keep_default_na=False, delimiter="\t")
    df_tag_spec["tag_len"] = df_tag_spec.apply(lambda x: len(str(x["tag"])), axis=1)
    df_tag_spec_sorted = df_tag_spec.sort_values(by=["tag_len"], ascending=False)
    #drop duplicate tags by position
    df_tag_spec_sortedNR=df_tag_spec_sorted.drop_duplicates(subset=["spectrum","tag_index"], keep="first")    
    spec_tag_df_cnt=df_tag_spec_sortedNR.groupby('spectrum').apply(my_agg).reset_index()
    
    spec_tag_df_cnt.to_csv("spectrum_unique_tag_table.txt", sep="\t", index=None)
    spectrum_tag_total_dict = dict(zip(spec_tag_df_cnt.spectrum,spec_tag_df_cnt["Total_Tag#"]))

    return spectrum_tag_total_dict

#number of matchd tags len(spectrum_tag_dict[key])
#key = spectrum+"_"+modifiedpeptide

spectrum_tag_total_dict = total_tags_compute(out_file)

#reorder_psms(dfMz_valid)
#function writePepXml in pepxml_Generation_Reorder file is changed to rank 1 and 2

out_tag_pep = "{}/{}_tag_matched.1.pep.xml".format(output_directory,fraction_name)
makePepXML(dfMz_valid, pep_xml,out_tag_pep, spectrum_tag_total_dict, logFile)



mv_cmd = "mv {} {}".format(out_tag_pep,pep_xml)
os.system(mv_cmd)

#weightedEvalue = combination of zscore(Escore = -log10(expect) and tagScore)
df = qc_target_decoy(dfMz_valid,logFile, "search_hit", "all", "num_matched_tags")
write_log (logFile,"  The num_matched_tags qc numbers are shown below")
qc_figures(df, "num_matched_tags", "Total_Tag_matched")


df = qc_target_decoy(dfMz_valid,logFile, "search_hit", "all", "expect")
write_log (logFile,"  The comet expect value qc numbers are shown below")
qc_figures(df, "-log10(eval)", "expect_minusLog10")


df = qc_target_decoy(dfMz_valid,logFile, "search_hit", "all", "xcorr")
write_log (logFile,"  The comet xcorr value qc numbers are shown below")
qc_figures(df, "xcorr", "xcorr")

write_log (logFile,"    Total Targets and Decoys ",df.Type.value_counts())



cmd_cp = "cp {}/{} .".format(inpath,params_file)
os.system(cmd_cp)


