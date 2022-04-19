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
import subprocess


from os.path import dirname

# path for the script
source_path_script = sys.argv[0]

#get the path of script
source_path = dirname(source_path_script)



def msg(name=None):
    return '''\n\n\npython ptm_pipeline.py jump ptm_pipeline.params file.mzXML/s\n\n\n '''

parser = argparse.ArgumentParser(description="Stagewise search using PTM Pipeline", prog="ptm_pipeline",usage=msg())
parser.add_argument("ptm_pipeline_parameterfile", help="ptm pipeline parameter file")
parser.add_argument("mzXML",help="single or list of mzXML files",nargs='+')
parser.add_argument("--queue",dest="queue", default="queue=standard")
parser.add_argument("--mem",dest="mem", default="mem=8000")

args = parser.parse_args()
mzXMLs_all = args.mzXML
queue = args.queue.split("=")[-1]
mem = args.mem.split("=")[-1]

mzXMLs = []
for mzF in mzXMLs_all:
    mzPath = os.getcwd()+"/{}".format(mzF)
    mzXMLs.append(mzPath)

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

jump_d_program = "perl {}/JUMPd/bin/_jump_d_v2.pl".format(source_path)
jump_qc_program = "perl {}/qc/spectrumQC.pl".format(source_path)
jump_q_program = "perl {}/JUMPq/jump_q.pl".format(source_path)
jump_q_params = "{}/parameterFiles/jump_qc_HH_tmt10_human.params".format(source_path)

#tag program
jump_tag_program = "python {}/JUMPtagmatch/main.py".format(source_path)
###### Store ptm_pipeline parameter file as the dictionary #####

allParamLines, allComments, allParamsDict  = storeParamsFile(params_file)


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
comet = "{}/comet".format(source_path)


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

#check for pitfile if pitfile == 0 it generates pitfile using jump -d
#jump -d use for generation of pit file and 

database_name = allParamsDict["database_name"]
pitfile = allParamsDict["pitfile"]

try:
    cmd1 = "rm -r database"
except:
    print ("  No previous pit file generated")

makedirectory("database")


if pitfile == "0": #need to apply jump -d to generate protein inference table (pit) file
    write_log(logFile,"jump -d program activated")
    write_log(logFile,"  The user did not specify the pit file so running jump -d to generate pitfile")
    #gathering parameters from ptm pipeline params file
    
    os.chdir("database")

    jump_d_params = "{}/parameterFiles/jump_d.params".format(source_path)

    input_database1 = allParamsDict["database_name"]
    output_prefix = allParamsDict["output_prefix"]
    include_contaminants = allParamsDict["include_contaminants"]
    decoy_generation = allParamsDict["decoy_generation"]
    decoy_generation_method = allParamsDict["decoy_generation_method"]
   
    generate_jump_d_param_file(jump_d_params, input_database1, output_prefix, include_contaminants, decoy_generation, decoy_generation_method)
    write_log(logFile,"  jump -d parameter file generated in the fly based on user parameters") 
    write_log(logFile,"  parameters used: \n  database_name\n  output_prefix\n  include_contaminants\n  decoy_generation\n  decoy_generation_method")

    cmd = "{} {}".format(jump_d_program,"jump_d_fly.params")
    os.system(cmd)


    pitfile = curr_path+"/{}/database/{}_ptm_pipeline_comet.pit".format(result_folder,output_prefix)
    database_name = curr_path+"/{}/database/{}_ptm_pipeline_comet.fasta".format(result_folder,output_prefix)

    allParamsDict["pitfile"] = pitfile
    allParamsDict["database_name"] = database_name

    write_log(logFile,"  New database to be used for Stage_0 searches {}".format(pitfile)) 
    write_log(logFile,"  New pitfile to be used for Stage_0 searches {}".format(database_name))



    os.chdir(curr_path+"/{}".format(result_folder))
    
else:
    write_log(logFile,"  The user supplied pitfile information {} so jump -d is skipped".format(pitfile))

#run comet params here to get the empty parameter file
subprocess.call([comet,"-p"])
write_log(logFile,"  The pipeline created masters comet parameter file -- comet.params.new and stored in the {} folder\n. This master file will be updated based on user input for different stages".format(os.getcwd())) 


comet_params = "comet.params.new"

#check mzXML files present or absent
if len(mzXMLs) == 0:
    sys.exit("FATAL: The program expects mzXML file in your current path. Please make sure you have mzXML in the current path")
    


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



stage0_folder = curr_path+"/{}/Stage_0".format(result_folder)
preprocess = curr_path+"/{}/Preprocess".format(result_folder)
makedirectory(preprocess)


if idtxt_file == "0":
    print ("    Please provide JUMP -f ID.txt file path with desired FDR = 1% (peptide level)")
    sys.exit()
    
else:
    write_log (logFile,"\n\nThis part of logFile consists of Stage_0 (whole unmodified search information")
    # stage0_folder = preprocess
    write_log (logFile,"  It seems like stage_0 whole proteome results is provided. The search results are in {} folders. Skipping Stage_0 automatic comet search. Thus no preprocessing required".format(preprocess))
    
    
    if allParamsDict["whole_proteome_path"] != "0":
        dtasFolder = allParamsDict["whole_proteome_path"]
        write_log (logFile,"  whole proteome path {} is provided so the program will not infer path from ID.txt file provided".format(dtasFolder))

    else:
        dtasFolder = find_stage0_folders_jump(idtxt_file) #this will mimic preprocess = curr_path+"/Pipeline_Results/Preprocess"
        write_log (logFile,"  whole proteome path is not provided so the program will infer path from ID.txt file provided. The infered path is {}".format(dtasFolder))


    #this means there is no modificaitions
    write_log (logFile,"  Since this is whole proteome there is no modified peptides. The program will now generate the ms2 files removing accepted psms from {}".format(idtxt_file))
    prev_mod = "0"
    
    all_dtas = []
    all_fractions_search = []
    #use mzXML file to get the desired fractions dtas
    for mz_file in mzXMLs:
        basefile = os.path.basename(mz_file)
        basefile_noext = basefile.split(".")[0]

        #add all fractions to the list. This is required for qc later if someone wants to work on the subset of mzXML files
        all_fractions_search.append(basefile_noext)
        # dta_file = glob.glob(dtasFolder+"/{}/{}.*dtas".format(basefile_noext,basefile_noext))[0]
        # all_dtas.append(dta_file)
    
    
    
    ####For QC spectrum analysis ####
    spectrum_qc_results_path = allParamsDict["spectrum_qc_results_path"]

    if spectrum_qc_results_path != "0":
        dtasFolder = spectrum_qc_results_path

    if (spectrum_qc_results_path == "0") and (enable_spectrum_QC == "1"):
    
        write_log (logFile,"\n\nThis part of logFile has spectrum QC information")
        write_log (logFile,"  User enabled spectrum QC so the pipeline will use spectrum_QC program. For more details of this program https://pubmed.ncbi.nlm.nih.gov/27225868/")
        write_log (logFile,"  Following parameters are internally defined for high confidence filtering")
        write_log (logFile,"  FDR=0\n  unique_protein_or_peptide=peptide\n  one_hit_wonders_removal=0\n  min_XCorr=15 (JUMP -search specific)\n  search_engine=jump\n  pitfile={}\n  database={}".format(pitfile, database_name))
        write_log (logFile,"  QC program part is activated. Here, the program generates FDR = 0 (high confidence) filter result at the peptide level.\nAlso, the program uses ID.txt file given by the user as the alternative filtering results.")
        
        os.chdir(stage0_folder)
        
        allParamsDict_f["FDR"] = "0"
        allParamsDict_f["unique_protein_or_peptide"] = "peptide"
        allParamsDict_f["one_hit_wonders_removal"] = "0"
        allParamsDict_f["min_XCorr"] = "15"
        allParamsDict_f["search_engine"] = "jump"
        allParamsDict_f["pit_file"] = pitfile
        allParamsDict_f["database"] = database_name


        #prepare_jump_f_input(folder, fdr, pitfile, stage, search_engine, min_XCorr)

        print ("    Dtas folder identified by program is {}".format(dtasFolder))
        out_params = prepare_jump_f_paramsFile(dtasFolder, allParamsDict_f, os.path.basename(stage0_folder), allParamsDict_f["FDR"])

        # out_params,list_pubs, fdr_add, pitfile_str,database_str, search_engine_str, min_XCorr_str = prepare_jump_f_input(dtasFolder, fdr, pitfile, database_name,os.path.basename(stage0_folder),"jump","10")#search engine and min
        # print ("    database  {}".format(database_str))
        # print ("    pit  {}".format(pitfile_str))
        # from_sample_to_batch_params(jump_f_params, out_params, list_pubs, fdr_add, pitfile_str,database_str, search_engine_str, min_XCorr_str)

        # cmd = "jump -f "+out_params
        cmd = "{} {}".format(jump_f_program,out_params)
        os.system(cmd)
        jump_f_id_txt_0 = "{}/sum_Stage_0_FDR_0/ID.txt".format(stage0_folder)

        #this is required to make the parameter file
        #confident_IDtxt, accepted_IDtxt
        qc_id_txt_list = [jump_f_id_txt_0,idtxt_file]

        #generate QC parameters
        write_log (logFile,"  The required high confidence FDR = 0 idtxt file and user provided idtxt file is used for QC program")
        generate_QC_param_file(dtasFolder, "jump_qc.params", qc_id_txt_list, PSM_recoveray_rate, all_fractions_search)
        write_log (logFile,"  QC parameter file jump_qc.params is generated in the fly and stored in {} folder".format(os.getcwd()))
        cmd = "{} {}".format(jump_qc_program,"jump_qc.params")
        os.system(cmd)

    for frac in all_fractions_search:
        if (spectrum_qc_results_path == "0") and (enable_spectrum_QC == "1"):
            dta_file = glob.glob("qc_accepted_PSMs/{}/{}*dtas".format(frac,frac))[0]
            all_dtas.append(dta_file)

        else:
            dta_file = glob.glob(dtasFolder+"/{}/{}*dtas".format(frac,frac))[0]
            all_dtas.append(dta_file)

    write_log (logFile,"  QC program completed and now using filtered outfiles to generate high quality spectrum ms2 files")
    for dtas in all_dtas:
        basefile = os.path.basename(dtas)
        basefile_noext = basefile.split(".")[0]
        makedirectory(preprocess+"/"+basefile_noext) 
        new_ms2 = preprocess+"/"+basefile_noext+"/"+basefile_noext+".ms2"
        write_log(logFile,"    Generating {} file".format(new_ms2))
        dta_to_ms2(dtas, new_ms2)
        write_log (logFile,"  The high quality ms2 file {} is written".format(new_ms2))



# #### MAKE MS2 REMOVING ACCEPTED PSMS ###########
os.chdir(stage0_folder)

write_log (logFile,"\n\nThis section of logFile has filtering information for Stage_0")

fdr = allParamsDict["FDR"]
jump_f_id_txt = idtxt_file

accepted_protein_list, scanChargeDf = idtxt_exclude_acceptedPSMs_acceptedProtein(jump_f_id_txt)


additional_qc_file = allParamsDict["additional_qc_file"]

if additional_qc_file != "0":
    qc_scan_df = additional_qc_parse(additional_qc_file, logFile)
    frames = [qc_scan_df, scanChargeDf]
    scanChargeDf = pd.concat(frames)

#generate a custom database
#make a folder called custom inside the database

if allParamsDict["custom_database"] == "1":
    write_log (logFile,"\n\nThis section of logFile has custom DB generation information")
    custom_db_folder = curr_path+"/{}/database/customDB".format(result_folder)
    makedirectory(custom_db_folder)

    write_log (logFile,"  The custom database based on filtered ID from Stage_0 is being generated")
    customDB = filtered_IDs_to_fasta(database_name,custom_db_folder, accepted_protein_list)
    database_name = customDB
    allParamsDict["database_name"] = database_name
    write_log (logFile,"  New custom DB is generated and stored as {} ".format(database_name))


#user is expected to have total ptm stages parameter same as the number of stage they want to search. This part checks and gives warnign if the stage number does not match
#this is not FATAl so the program goes on
if check_stage < int(allParamsDict["total_ptm_stages"]):
     write_log (logFile,"   \n\nWARNING !!! Total stages of PTMs defined by user does not match with the Stagewise parameters.") 
     write_log (logFile,"   This does not affect the search but if you meant to search {} stages; currently you only specified {} stages parameters".format(allParamsDict["total_ptm_stages"],check_stage)) 


os.chdir(curr_path+"/{}".format(result_folder))

write_log (logFile,"\n\nThis section of logFile has multiple stages PTM searches information")


ptms_list = []
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

            #gives the list of ptms that are being searched. 
            working_ptms = all_stages_parameters_dict["ptm_stage_{}".format(x)]
            ptms_list.append(working_ptms)

            write_log (logFile," This stage will search {} ptms".format(working_ptms))

            #updates check_stage with each stages
            check_stage+=1

            #this updates the stagewise parameters and updates the comet dictionary
            stagewise_var_mods_update(x, all_stages_parameters_dict, cometParamsDict, variable_mod01_mod_number)

            #this make stagewise parameter file for comet ptm searching     
            makeCometParametersStageWise(ptm_params, cometParamLines, cometComments, cometParamsDict)

            write_log (logFile," PTM parameter file for {} ptms is generated and stored as {} ".format(working_ptms, ptm_params))


write_log (logFile,"  List of all PTMs to be searched")
write_log (logFile,"    {}".format("\n  ".join(ptms_list)))

print ("Total stages to search = \n{}".format(stages_folder))

for folders in stages_folder:
    if "Stage_0" not in folders:
        stage_wise_search(folders, mzXMLs, comet, preprocess, scanChargeDf)
        # stage_wise_search_server(folders, mzXMLs, comet,logFile, preprocess, scanChargeDf)
        # wait_Function(1, folders, mzXMLs)
        rename_pep_xml(mzXMLs, folders)

        run_tag_program(mzXMLs, folders, jump_tag_program, dtasFolder)

        # wait_function_after_tag_match(mzXMLs, folders, "tag_qc.params")


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

allParamLines_q, allComments_q, allParamsDict_q  = storeParamsFile(jump_q_params)

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

#merge idtxt_file as Stage_0 jump -f results 

write_log (logFile,"  Concatenating jump -f ID.txt")
head1,df_idtxt_all = combine_jump_f_stages_modPeptides(all_stage_idtxts, idtxt_file)


combine_filter_result = allParamsDict["combine_filter_result"]


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


write_log (logFile,"  Concatenating jump -f id_uni_pep.txt")

#use stage0_idtxt to find id_uni_pep file and merge from stage0
df_uni_pep = merge_publication_jump_f_tables(all_stage_id_uni_pep, merge_directory, idtxt_file, keep_spectrum)
write_log (logFile,"  Concatenating jump -f id_all_pep.txt")
df_all_pep = merge_publication_jump_f_tables(all_stage_id_all_pep, merge_directory, idtxt_file, keep_spectrum)

#Options for QC after concatenation of jump -f stagewise results: 1: Xcorr based priority (best xcorr scan retains) [DEFAULT], 2: Stage priority (Stage 1 wins over Stage2), 3: Keep all result (no removal)- Chances few scans could have different peptides but identical quantification), 4: Remove all redundant scans

write_log (logFile,"  Concatenation of ID.txt files, id_uni_pep.txt files and id_all_pep.txt files successful. The results are stored in {}".format(os.getcwd()+"/merge_all_stages"))



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


###### Combination of jump -f and perform jump -q on the concatenated results #####
allParamLines_q, allComments_q, allParamsDict_q  = storeParamsFile(jump_q_params)

os.chdir(merge_directory)

write_log (logFile,"\n\nThis section of logFile has jump -q quantification results on merged jump -f results and stored in {}".format(os.getcwd()))

#merge stage_0 idtxt file for proper jump -q loading bias results



allParamsDict_q["idtxt"] = curr_path+"/{}/{}/ID.txt".format(result_folder,merge_directory)
out_params = prepare_jump_q_paramsFile(allParamsDict_q, "merged")
cmd = "{} {}".format(jump_q_program,out_params)

os.system(cmd)

write_log(logFile,"All ptm stages are successfully searched, filtered, merged and quantified.")



