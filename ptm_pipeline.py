import sys, os
import pandas as pd
import glob

import argparse
import numpy as np
import re
from logFunctions import *
from parameterFileFunctions import *
from pep_xml_parser import *
from stagewise_comet_search import *
from format_pep_xml import *
from idtxt_merge_filter import *
from customDB_gen import *
from program_signature import program_display
from postQuantificationFunctions import *
from id_all_uni_from_id import *
import subprocess
import platform

from os.path import dirname

# path for the script
source_path_script = sys.argv[0]

#get the path of script
source_path = dirname(source_path_script)



def msg(name=None):
    return '''\n\n\npython ptm_pipeline.py jump ptm_pipeline.params file.mzXML/s\n\n\n '''

parser = argparse.ArgumentParser(description="Stagewise search using PTM Pipeline", prog="ptm_pipeline",usage=msg())
parser.add_argument("ptm_pipeline_parameterfile", help="ptm pipeline parameter file")
parser.add_argument("--queue",dest="queue", default="queue=standard")
parser.add_argument("--mem",dest="mem", default="mem=8000")

args = parser.parse_args()
queue = args.queue.split("=")[-1]
mem = args.mem.split("=")[-1]



curr_path = os.getcwd()

params_file = curr_path+"/"+args.ptm_pipeline_parameterfile
print ("Parameter file is {}".format(params_file))
#params_file = "/home/spoudel1/PanPTM_Paper_2021/Pipeline/parameterFiles/ptm_pipeline_v001_YL.params"

#Change the directory to working directory

#we will run the program in the directory that contains the mzxml files

# working_path = "/home/spoudel1/PanPTM_Paper_2021/Figure1/Test_PTM_Pipeline"
# os.chdir(working_path)

#mzXMLs = glob.glob(curr_path+"/*.mzXML")


##### All required programs and parameter files for pipeline #####
#### 1. JUMP -f program and its parameter file
#### 2. Deisotope program 
#### 3. JUMP -d program
#### 4. JUMP -spectrum_QC program
#### 5. JUMP -q program and its parameter file


jump_f_program = "perl {}/JUMPf/jump_f.pl".format(source_path)
#tag based filteration
#jump_f_program = "perl /home/yli4/development/git/Yuxin_jump_git/JUMPf/jump_f.pl"
# jump_f_program = "perl /hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.004/JUMPf/bin/_jump_f.pl"
#deisotope_program = "python {}/Deisotope/DeisotopeMS1.py".format(source_path)
jump_f_params = "{}/parameterFiles/jump_fc_HH_tmt10_human.params".format(source_path)
jump_q_program = "python {}/JUMPq/jump_q.py".format(source_path)
jump_q_params = "{}/parameterFiles/jump_qc_HH_tmt10_human.params".format(source_path)

#tag program
jump_tag_program = "python {}/JUMPtagmatch/main.py".format(source_path)
###### Store ptm_pipeline parameter file as the dictionary #####

allParamLines, allComments, allParamsDict  = storeParamsFile(params_file)

ms2_path = allParamsDict["ms2_input_path"]
tags_input_path = allParamsDict["tags_input_path"]

mzXMLs = glob.glob(ms2_path+"/*.ms2")


idtxt_file = allParamsDict["idtxt_file"]

#some mandatory parameters will be overwritten by the program
#for example generation of search .txt files and pep.xml files are expected so we will give those parameters as internal parameters
allParamsDict["output_txtfile"] = "1"
allParamsDict["output_pepxmlfile"] = "1" 
allParamsDict["output_suffix"] = ""



#This is for QC spectrum
enable_spectrum_QC = allParamsDict["enable_spectrum_QC"]
PSM_recoveray_rate = allParamsDict["PSM_recoveray_rate"]

# Comet 2021 version binary file
if platform.system() == "Linux":
    comet = "{}/comet_linux_2021".format(source_path)

if platform.system() == "Darwin":
    comet = "{}/comet_mac_2021".format(source_path)


database_name = allParamsDict["database_name"]
pitfile = allParamsDict["pitfile"]

###extract PTM pipeline dictionary from master dictionary
all_stages_parameters_dict = {}


##only stage related parameters
for all_parameters in allParamsDict.keys():
    if "ptm_stage" in all_parameters or "max_mods_per_peptide_stage" in all_parameters: #ptm pipeline parameter found
        all_stages_parameters_dict[all_parameters] = allParamsDict[all_parameters]
        

#this sets the current directory as the working directory. So any new files generated will be saved in this folder
#extract suffix to be added to the output folder 
folder_suffix = allParamsDict["ptm_output_folder_suffix"] 
result_folder = "Pipeline_Results_{}".format(folder_suffix)
makedirectory(result_folder)


#copies ptm_parameter file to results folder
cpFile(params_file, result_folder)


os.chdir(curr_path+"/{}".format(result_folder))
logFile = curr_path+"/{}/ptm_pipeline.log".format(result_folder)
#rmFile("ptm_pipeline.log")
try:
    rmFile(logFile)
except:
    write_log(logFile, "No previous logFile present for PTM. So a new logFile is generated")


#### Launch the program ####
display_Text = program_display()
write_log(logFile,display_Text)

write_log (logFile," ********* PTM Pipeline Initiated **************")


#run comet params here to get the empty parameter file
subprocess.call([comet,"-p"])
write_log(logFile,"  The pipeline created masters comet parameter file -- comet.params.new and stored in the {} folder\n. This master file will be updated based on user input for different stages".format(os.getcwd())) 


comet_params = "comet.params.new"

#check mzXML files present or absent
if len(mzXMLs) == 0:
    sys.exit("FATAL: The program expects ms2/mzXML file in your current path. Please make sure you have ms2 in the current path")
    


#use this part to generate new parameter files based on the ptm pipeline parameter file and populate the empty comet paramter file
#extract the methionine oxidation total mods per peptide if described by users
#this will be added to total modification per stage to design final modification number per peptide
variable_mod01_mod_number = int(allParamsDict["variable_mod01"].split(" ")[3])
write_log (logFile,"\n\nThis section of log file has the stagewise information")
write_log (logFile,"   A total of {} PTM stages are given by the user for searching".format(allParamsDict["total_ptm_stages"]))

#this counts the number of valid stages
check_stage = 0
#total number of ptms stages = allParamsDict["total_ptm_stages"])

stages_folder = []
#make stage 0 folders
makedirectory("Stage_0")
stages_folder.append(os.getcwd()+"/Stage_0")

write_log (logFile," \n\n****** Generating Parameter files *********\n")
write_log (logFile,"  Working for stage 0")
write_log (logFile,"  Stage_0 results are provided. So, the program will use unmodified jump -s for Stage_0")



#JUMP -f on the Stage_0 search results

allParamLines_f, allComments_ff, allParamsDict_f  = storeParamsFile(jump_f_params)




os.chdir(curr_path+"/{}".format(result_folder))

write_log (logFile,"\n\nThis section of logFile has multiple stages PTM searches information")

stage_ptm_dict = {}
ptms_list = []
dynamic_mods_aa_dict = {}
ptm_aa_dict = {}
for x in range(1,int(allParamsDict["total_ptm_stages"])+1):

    if "ptm_stage_{}_amino_acids".format(x) in allParamsDict.keys():
        if allParamsDict["ptm_stage_{}_amino_acids".format(x)] != "X":

            #every stage directory is made
            makedirectory("Stage_"+str(x))

            stages_folder.append(os.getcwd()+"/Stage_"+str(x))

           
            #store the empty parameter file from comet
            #this should be updated for every stages
            cometParamLines, cometComments, cometParamsDict = storeParamsFile(comet_params)
            #update the new comet parameter dictionary with user defined all ptm pipeline parameter dictionary
            #only same keys are updated
            update_comet_dictionary(cometParamsDict, allParamsDict)


            #comet parameter file is formed inside the same stage using the stagewise parameters
            ptm_params = "Stage_{}/stage_{}_comet.params".format(x,x)

            write_log (logFile," \n\n****** Generating Parameter files *********\n")
            write_log (logFile," Working for stage {}".format(x))

            
            ptms_info_params(all_stages_parameters_dict, x, stage_ptm_dict, dynamic_mods_aa_dict,ptm_aa_dict)

            working_ptms = stage_ptm_dict["Stage_"+str(x)]
            write_log (logFile," This stage will search {} ptms".format(working_ptms))

            #updates check_stage with each stages
            check_stage+=1

            #this updates the stagewise parameters and updates the comet dictionary
            stagewise_var_mods_update(x, all_stages_parameters_dict, cometParamsDict, variable_mod01_mod_number)

            #this make stagewise parameter file for comet ptm searching     
            makeCometParametersStageWise(ptm_params, cometParamLines, cometComments, cometParamsDict)

            write_log (logFile," PTM parameter file for {} ptms is generated and stored as {} ".format(working_ptms, ptm_params))



# print (stage_ptm_dict)
# print (dynamic_mods_aa_dict)
# print(ptm_aa_dict)

#user is expected to have total ptm stages parameter same as the number of stage they want to search. This part checks and gives warnign if the stage number does not match
#this is not FATAl so the program goes on
if check_stage < int(allParamsDict["total_ptm_stages"]):
     write_log (logFile,"   \n\nWARNING !!! Total stages of PTMs defined by user does not match with the Stagewise parameters.") 
     write_log (logFile,"   This does not affect the search but if you meant to search {} stages; currently you only specified {} stages parameters".format(allParamsDict["total_ptm_stages"],check_stage)) 


write_log (logFile,"  List of all PTMs to be searched")
write_log (logFile,"    {}".format("\n  ".join(ptms_list)))

#print ("Total stages to search = \n{}".format(stages_folder))

scanChargeDf=None
cluster = allParamsDict["cluster"]

for folders in stages_folder:
   if "Stage_0" not in folders:
        stage_wise_search(folders, mzXMLs, comet, logFile, cluster, scanChargeDf)
        rename_pep_xml(mzXMLs, folders)
        if tags_input_path == "0":
            write_log(logFile,"  Tags file are missing so the JUMPptm will be performed without tag based filtering")
        else:
            run_tag_program(mzXMLs, folders, jump_tag_program, tags_input_path, cluster)
   

# #### JUMP -f filter on each stages #####

write_log (logFile,"\n\nThis section of logFile has multiple stages jump -f information")

fdr = allParamsDict["FDR"]

for folders in stages_folder:
    if os.path.basename(folders) != "Stage_0":
        os.chdir(folders)

            # allParamsDict_f["FDR"] = "0"
        # allParamsDict_f["unique_protein_or_peptide"] = "peptide"
        # allParamsDict_f["one_hit_wonders_removal"] = "0"
        allParamsDict_f["min_XCorr"] = "1"
        allParamsDict_f["search_engine"] = "comet"
        allParamsDict_f["pit_file"] = pitfile
        allParamsDict_f["database"] = database_name

        # if idtxt_file != "0":
        #     out_params = prepare_jump_f_paramsFile_tag(folders, allParamsDict_f, os.path.basename(folders), fdr)
        # else:
        out_params = prepare_jump_f_paramsFile(folders, allParamsDict_f, os.path.basename(folders), fdr)

#         cmd = "jump -f "+out_params
        write_log (logFile,"  Running jump -f for {}".format(os.getcwd()))
        cmd = "{} {}".format(jump_f_program,out_params)

        os.system(cmd)

###### Combination of jump -f and perform jump -q on the concatenated results #####

write_log (logFile,"\n\nThis section of logFile has concatenation of jump -f results")



all_stage_idtxts = []
all_stage_id_uni_pep = []
all_stage_id_all_pep = []

for folders in stages_folder:
    if os.path.basename(folders) != "Stage_0":
        idtxt_Stage = glob.glob(folders+"/sum_Stage_*/ID.txt")[0] 
        uniq_pep_file = glob.glob(folders+"/sum_Stage_*/publications/id_uni_pep.txt")[0] 
        all_pep_file = glob.glob(folders+"/sum_Stage_*/publications/id_all_pep.txt")[0] 
        
        all_stage_idtxts.append(idtxt_Stage)
        all_stage_id_uni_pep.append(uniq_pep_file)
        all_stage_id_all_pep.append(all_pep_file)


os.chdir(curr_path+"/{}".format(result_folder))
#make a out_folder for merged jump _f results
merge_directory = "merge_and_consolidation"
makedirectory(merge_directory)
makedirectory(merge_directory+"/publications")



#merge idtxt_file as Stage_0 jump -f results 

# print ("      Stage 0 idtxt = {}\n\n".format(stage0_idtxt))
write_log (logFile,"  Concatenating jump -f ID.txt")
head1,df_idtxt_all = combine_jump_f_stages_modPeptides(all_stage_idtxts, idtxt_file)


combine_filter_result = allParamsDict["combine_filter_result"]


keep_spectrum = []

if combine_filter_result == "1":
    #we get the list of spectra (stage_spectrum) to be retained based on 
    write_log (logFile, "  You selected keep one redundant spectra that has highest XCorr if occured in multiple stages. So, the output will be merged idtxt minus all duplicate spectra (retaining one based on XCorr)")
    keep_spectrum=get_non_redundant_spectra_xcorr(df_idtxt_all)

if combine_filter_result == "2":
    write_log (logFile, "  You selected keep one redundant spectra that belongs to the earliest stage. So, the output will be merged idtxt minus all duplicate spectra (retaining one based on stage)")
    keep_spectrum=get_non_redundant_spectra_stagewise(df_idtxt_all)

    
if combine_filter_result == "4":
    write_log (logFile, "  You selected remove all redundant spectra from different stages if any. So, the output will be merged idtxt minus all redundant spectra")
    keep_spectrum=get_non_redundant_spectra_remove_all_duplicates(df_idtxt_all)

if combine_filter_result == "3":
    write_log (logFile, "  You selected not to remove redundant spectra from different stages if any. So, the output will be merged idtxt")

df_idtxt_all = combine_filter_result_consolidate(df_idtxt_all,keep_spectrum)

write_log (logFile,"  filtering option combine_filter_result = {} applied to prepare ID.txt for jump -q".format(combine_filter_result))

#now process this merged idtxt file to generate the ID.txt file in the merged folder
#save the idtxt file without the header of datbase
df_idtxt_all.to_csv(merge_directory+"/ID_no_header.txt",sep=";", index=None)


ori_clean1 = merge_directory+"/ID_no_header.txt"

lines1=[]
g1=open(ori_clean1,"r")
lines1 = g1.readlines()
lines1.insert(0,head1.strip())

with open(merge_directory+"/ID.txt","w") as final1:
    for value1 in lines1:
        final1.write(value1.strip()+"\n")


generate_id_uni_all_txt_files(pitfile, df_idtxt_all, merge_directory)

###### Combination of jump -f and perform jump -q on the concatenated results #####

allParamLines_q, allComments_q, allParamsDict_q  = storeParamsFile(jump_q_params)
# replace impurity matrix from github resources impurity_matrix =
allParamsDict_q["impurity_matrix"] = "{}/resources/{}".format(source_path,allParamsDict_q["impurity_matrix"])
print ("impurity matrix path is ",allParamsDict_q["impurity_matrix"])

os.chdir(merge_directory)

write_log (logFile,"\n\nThis section of logFile has jump -q quantification results on merged jump -f results and stored in {}".format(os.getcwd()))

#merge stage_0 idtxt file for proper jump -q loading bias results



allParamsDict_q["idtxt"] = curr_path+"/{}/{}/ID.txt".format(result_folder,merge_directory)
out_params = prepare_jump_q_paramsFile(allParamsDict_q, "merged")
cmd = "{} {}".format(jump_q_program,out_params)

os.system(cmd)

write_log(logFile,"All ptm stages are successfully searched, filtered, merged and quantified.")

#get the final quantification matrix
# print (result_folder,"\n",merge_directory)
in1 = glob.glob("{}/{}/{}/quan_*/publications/id_uni_pep_quan.txt".format(curr_path,result_folder,merge_directory))[0]
# in1 = "/research_jude/rgs01_jude/groups/penggrp/projects/Alzheimer_BannerSunInstitute/penggrp/proteomics/batch0_pooledSamples/PanPTM_Paper_2021/github_test_6_small_file/Pipeline_Results_jumpptm_simplified/merge_and_consolidation/quan_HH_tmt10_human_comet/publications/id_uni_pep_quan.txt"



write_log(logFile,"Generating final pan ptm excel table. This will take time depending on the number of quantified peptides")
#final quantification matrix from jump -q

df_batch0_q = pd.read_csv(in1, delimiter="\t", low_memory=False)

#stage_ptm_dict = {'Stage_1': 'phosphorylation', 'Stage_2': 'Deamidation'}


channels_cols = []
for x in df_batch0_q.columns:
    if "sig" in x:
        channels_cols.append(x)

keep_cols = ['Peptides', 'Peptides_original','Protein Group#', 'Protein Accession #',
       'Protein Description', 'GN', 'PSM#', 'Run#', 'Scan#','m/z', 'z',
       'XCorr',  'ptm_stage'] + channels_cols

df_batch0_q["PTM_types"] = df_batch0_q.ptm_stage.map(stage_ptm_dict)
#result_folder = dirname(dirname(dirname(dirname(in1))))



all_pep_xmls = glob.glob("{}/{}/Stage*/*/*.pep.xml".format(curr_path,result_folder))
stage_pep_dict = {}
for peps in all_pep_xmls:
    stage = os.path.basename(dirname(dirname(peps)))
    if stage not in stage_pep_dict.keys():
        stage_pep_dict[stage] = peps


#dynamic_mods_aa_dict={79.966331: ['S', 'T', 'Y'], 0.984016: ['N', 'Q']}
#ptm_aa_dict = {'Stage_1': ['S', 'T', 'Y'], 'Stage_2': ['N', 'Q']}

replace_symbol_with_val(df_batch0_q, stage_pep_dict, dynamic_mods_aa_dict)
make_spectrum_single_jump_f(df_batch0_q)

stage=os.path.basename(dirname(dirname(all_pep_xmls[0])))

df_0=extract_spectrum_eval_tag(all_pep_xmls[0], stage, tags_input_path)
for pep_file in all_pep_xmls[1:]:
    stage1=os.path.basename(dirname(dirname(pep_file)))
    df_1=extract_spectrum_eval_tag(pep_file,stage1,tags_input_path)
    frames = [df_0, df_1]
    df_0 = pd.concat(frames)


quick_row_iterate(df_0)

df_batch0_q["spectrum_stage_key"]=df_batch0_q["spectrum"]+"__"+df_batch0_q["ptm_stage"]
df_0["spectrum_stage_key"]=df_0["pseudo_spectrum"]+"__"+df_0["stage"]
df_batch0_q_eval_tag = df_batch0_q.merge(df_0, how="left", on="spectrum_stage_key")

if tags_input_path == "0":
    df_batch0_q_eval_tag["final_merge_key"] = df_batch0_q_eval_tag["spectrum__mod_pep"]+"__"+df_batch0_q_eval_tag["stage"]
    final_publ_cols = ['Peptides_with_ptm_mass','PSM#','expect','XCorr','Protein Accession #','GN','PTM_types']+channels_cols
    df_batch0_q_complete = df_batch0_q_eval_tag[final_publ_cols]
    rename_cols = {'Peptides_with_ptm_mass':'Peptides','PSM#':'# PSMs','expect':'Matching E values','Protein Accession #':'Protein ID','PTM_types':'PTM types'}


    df_batch0_q_final = df_batch0_q_complete.rename(columns=rename_cols)
else:

    unique_tag_files = glob.glob("{}/{}/Stage_*/*/Results*/spectrum_unique_tag_table.txt".format(curr_path,result_folder))
    tag_merge_df = merge_tag_files(unique_tag_files)

    select_best_tag(tag_merge_df)


    req_tag_columns = ["spectrum_stage_key","spectrum__mod_pep","stage","Top_scored_tag","Best_tag_scores"]
    tag_merge_df_table = tag_merge_df[req_tag_columns]

    #make final merge key
    df_batch0_q_eval_tag["final_merge_key"] = df_batch0_q_eval_tag["spectrum__mod_pep"]+"__"+df_batch0_q_eval_tag["stage"]
    tag_merge_df_table["final_merge_key"] = tag_merge_df_table["spectrum__mod_pep"]+"__"+tag_merge_df_table["stage"]

    #df_batch0_q_complete = df_batch0_q_eval_tag.merge(tag_merge_df_table, how="left", on = "spectrum_stage_key")
    df_batch0_q_complete = df_batch0_q_eval_tag.merge(tag_merge_df_table, how="left", on = "final_merge_key")


    final_publ_cols = ['Peptides_with_ptm_mass','PSM#','no of matched tags','Top_scored_tag', 'Best_tag_scores','expect','XCorr','Protein Accession #','GN','PTM_types']+channels_cols
    df_batch0_q_final = df_batch0_q_complete[final_publ_cols]


    rename_cols = {'Peptides_with_ptm_mass':'Peptides','PSM#':'# PSMs','no of matched tags':'# Matched tags',
                   'Top_scored_tag':'Top scored tag','Best_tag_scores':'Best tag scores','expect':'Matching E values','Protein Accession #':'Protein ID','PTM_types':'PTM types'}


    df_batch0_q_final = df_batch0_q_final.rename(columns=rename_cols)


df_batch0_q_final.drop_duplicates(inplace=True)
# make final table folder
makedirectory("results_table")
df_batch0_q_final.to_excel("results_table/Pan_PTM_Quan_Table.xlsx",index=None)


