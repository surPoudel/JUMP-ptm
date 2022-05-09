import configparser
import glob
config = configparser.ConfigParser()
from logFunctions import write_log, rmFile
import re
from os.path import dirname
import pandas as pd 
from stagewise_comet_search import return_skiprows

#we will not use this fucntion as the storeParamsFile works for reading the parameter file
def read_parameter_file(params_file):
    
    config.read(params_file)

    mzfile_path = config["ptm_pipeline"]["mzfile_path"].split("#")[0].strip()
    database = config["ptm_pipeline"]["database"].split("#")[0].strip()
    idtxt_file = config["ptm_pipeline"]["idtxt_file"].split("#")[0].strip()
    columns_names = config["ptm_pipeline"]["columns_names"].split("#")[0].strip()
    enable_spectrum_QC = config["ptm_pipeline"]["enable_spectrum_QC"].split("#")[0].strip()
    PSM_recoveray_rate = config["ptm_pipeline"]["PSM_recoveray_rate"].split("#")[0].strip()
    total_ptm_stages = int(config["ptm_pipeline"]["total_ptm_stages"].split("#")[0].strip())
    
    parameters_list=[mzfile_path,database,idtxt_file,columns_names,enable_spectrum_QC,PSM_recoveray_rate]
    
    all_stages_parameters_dict = {}
    for x in range(1,total_ptm_stages+1):
        ptm_stage = config["ptm_pipeline"]["ptm_stage_"+str(x)].split("#")[0].strip()
        ptm_stage_amino_acids = config["ptm_pipeline"]["ptm_stage_{}_amino_acids".format(x)].split("#")[0].strip()
        ptm_stage_dyn_mods = config["ptm_pipeline"]["ptm_stage_{}_dyn_mods".format(x)].split("#")[0].strip()
        max_mods_per_peptide_stage = config["ptm_pipeline"]["max_mods_per_peptide_stage_{}".format(x)].split("#")[0].strip()

        all_stages_parameters_dict["ptm_stage_{}".format(x)]=ptm_stage
        all_stages_parameters_dict["ptm_stage_{}_amino_acids".format(x)]=ptm_stage_amino_acids
        all_stages_parameters_dict["ptm_stage_{}_dyn_mods".format(x)]=ptm_stage_dyn_mods
        all_stages_parameters_dict["max_mods_per_peptide_stage_{}".format(x)]=max_mods_per_peptide_stage
        
    return all_stages_parameters_dict, parameters_list



def ptms_info_params(parameters_dict, stage, stage_ptm_dict, dynamic_mods_aa_dict, ptm_aa_dict):
    working_ptms = parameters_dict["ptm_stage_{}".format(stage)]
    dyn_mods = parameters_dict["ptm_stage_{}_dyn_mods".format(stage)]
    aa_residues = parameters_dict["ptm_stage_{}_amino_acids".format(stage)]
    
    working_ptms_list = working_ptms.split(";")
    dyn_mods_list = dyn_mods.split(";")
    aa_residues_list = aa_residues.split(";")

    for ind,mods in enumerate(dyn_mods_list):
        dynamic_mods_aa_dict[float(mods)] = list(aa_residues_list[ind])
    
    stage_ptm_dict["Stage_"+str(stage)] = ",".join(working_ptms_list)
    ptm_aa_dict["Stage_"+str(stage)] = dynamic_mods_aa_dict[float(mods)]


def storeParamsFile(paramFile):
    cometComments = {}
    dict1 = {}
    fileLines = []
    with open(paramFile,"r") as f:
        for line in f:
            if (line.startswith("#") == False) and ("=" in line.strip()):
                te_line = line.strip().split("=")
                key = te_line[0].strip()
                if "#" in te_line[1]:
                    te_line2 = te_line[1].split("#")
                    value1 = te_line2[0].strip()
                    valueComments = line.strip().split("#")[1]
                    cometComments[key] = "#"+valueComments
                else:
                    cometComments[key] = ""
                    value1 = te_line[1].strip()
                dict1[key]= value1
                fileLines.append(key)
            else:
                fileLines.append(line.strip())
    
    return fileLines,cometComments,dict1


def prepare_jump_f_paramsFile_tag(folder, allParamsDict, stage, fdr):
    out_params = "jump_fc_{}_FDR_{}.params".format(stage, "_".join(fdr.split(".")))#replace decimal with _ if any

    with open(out_params,"w") as f:
        batch_f_pub = glob.glob(folder+"/*/*_tag_matched.1.pep.xml")

        if len(batch_f_pub) == 0:
            print ("  Comet result is not there so checking for pepXML file")
            batch_f_pub = glob.glob(folder+"/*/*.pepXML")

        add_string="{}_FDR_{}:".format(stage, fdr)
        for pepfiles in batch_f_pub:
            # out_params = folder+"/jump_fc_stage_{}_FDR_.params".format(stage, fdr)
            
            in_string = pepfiles.split(".pep")[0]
            list_pubs = ""
            add_string+=in_string+"\n"
            list_pubs+=add_string

        f.write(list_pubs)
        
        for keys in allParamsDict.keys():
            f.write(keys+" = "+allParamsDict[keys]+"\n")
            
    return out_params


def prepare_jump_f_paramsFile(folder, allParamsDict, stage, fdr):
    out_params = "jump_fc_{}_FDR_{}.params".format(stage, "_".join(fdr.split(".")))#replace decimal with _ if any

    with open(out_params,"w") as f:
        batch_f_pub = glob.glob(folder+"/*/*.pep.xml")

        if len(batch_f_pub) == 0:
            print ("  Comet result is not there so checking for pepXML file")
            batch_f_pub = glob.glob(folder+"/*/*.pepXML")

        add_string="{}_FDR_{}:".format(stage, fdr)
        for pepfiles in batch_f_pub:
            # out_params = folder+"/jump_fc_stage_{}_FDR_.params".format(stage, fdr)
            
            in_string = pepfiles.split(".pep")[0]
            list_pubs = ""
            add_string+=in_string+"\n"
            list_pubs+=add_string

        f.write(list_pubs)
        
        for keys in allParamsDict.keys():
            f.write(keys+" = "+allParamsDict[keys]+"\n")
            
    return out_params

def prepare_jump_q_paramsFile(allParamsDict, stage):
    out_params = "jump_fq_{}.params".format(stage)
    with open(out_params,"w") as f:
        for keys in allParamsDict.keys():
            f.write(keys+" = "+allParamsDict[keys]+"\n")
    return out_params


def update_comet_dictionary(comet_params_dict, all_params_dict):
    for keys in all_params_dict.keys():
        if keys in comet_params_dict.keys():
            comet_params_dict[keys] = all_params_dict[keys]
            
            

            
def makeCometParametersStageWise(ptm_params, cometParamLines, cometComments, cometParamsDict):
    with open(ptm_params, "w") as paramC:
        for line in cometParamLines:
            if line.startswith("#"):
                paramC.write(line+"\n")
    #             print (line)      
            elif line == "":
                paramC.write(line+"\n")

            else:
                try:
                    paramC.write(line+" = "+cometParamsDict[line]+"\t"+cometComments[line]+"\n")
    #                 print (line," = ",cometParamsDict[line],"\t",cometComments[line])
                except:
                    paramC.write(line+"\n")


#ptm_stage_1 = phosphorylation
#ptm_stage_1_amino_acids = STY 
#ptm_stage_1_dyn_mods = 79.966331 
#max_mods_per_peptide_stage_1 = 3

def stagewise_var_mods_update(stage, all_stages_parameters_dict, cometParamsDict, variable_mod01_mod_number):
    #this dictionary is updated every stages 
    #generate a default stage 0 
    
    
    
    if stage == "0":
        cometParamsDict["max_variable_mods_in_peptide"] = "3"
        cometParamsDict["variable_mod02"] = "0.0 X 0 3 -1 0 0 0.0"
    else: 
        stage_dict = {}
        #only the parameters with stage_information are used to make stage_dict for every stage
        for stageKeys in all_stages_parameters_dict.keys():
            if "stage_{}".format(stage) in stageKeys:
                stage_dict[stageKeys] = all_stages_parameters_dict[stageKeys]

        #extracts the information for variable amino acid, modification mass and number of modificaition per peptides
        variable_aa_residues_list = stage_dict["ptm_stage_{}_amino_acids".format(stage)].split(";")
        variable_mods_mass_list =  stage_dict["ptm_stage_{}_dyn_mods".format(stage)].split(";")
        max_variable_mods_in_peptide_list = stage_dict["max_mods_per_peptide_stage_{}".format(stage)].split(";")


        #this step checks if user have properly provided the stagewise parameters; if not true the program exits with and error message
        if (len(variable_aa_residues_list)!=len(variable_mods_mass_list)) | (len(variable_aa_residues_list) !=len(max_variable_mods_in_peptide_list)):
            sys.exit("FATAL: Your ptm_stage_{}_amino_acids or ptm_stage_{}_dyn_mods  or max_mods_per_peptide_stage_{} are not even. Total entries separated by ; should be equal for these three parameters".format(stage,stage,stage))

        #max_variable_mods_in_peptide_per_stage is total stagewise modification number described by user + 2
        max_variable_mods_in_peptide_per_stage = sum(map(int, max_variable_mods_in_peptide_list))+variable_mod01_mod_number + 2
    #     print ("max_variable_mods_in_peptide_per_stage = ", max_variable_mods_in_peptide_per_stage)
        #updates maximum variable modifications per peptide per stage
        cometParamsDict["max_variable_mods_in_peptide"] = str(max_variable_mods_in_peptide_per_stage)


        #update comet params dictionary based on the stage information
        #initiate start value as 1 = oxidation of methionine
        #variable modification update starts from number 2
        #if user does not define methionine oxidation, the modificaion will be 0 for var mod 1
        for index, mass in enumerate(variable_mods_mass_list):
            aa = variable_aa_residues_list[index]
            max_mods = max_variable_mods_in_peptide_list[index]
            variable_mod_info = "{} {} {} {} {} {} {}".format(mass, aa, "0", max_mods, "-1", "0","0")

            #update comet dictionary
            cometParamsDict["variable_mod0{}".format(index+2)] = variable_mod_info

    return cometParamsDict


#this if for jump -f program parameter file generation automatically

def from_sample_to_batch_params(sample_params, output_params, replace_line1, replace_line2, replace_line3, replace_line4, replace_line5, replace_line6):
    with open(sample_params,"r") as f, open(output_params,"w") as g:
        for line in f:
            if "HH_tmt10_human_comet:" in line.strip():
                g.write(replace_line1+"\n")
                continue
            if "FDR =" in line.strip():
                g.write(replace_line2+"\n")
                continue
            if "pit_file = " in line.strip():
                g.write(replace_line3+"\n")
                continue

            if "database = " in line.strip():
                g.write(replace_line4+"\n")
                continue

            if "search_engine = " in line.strip():
                g.write(replace_line5+"\n")
                continue


            if "min_XCorr = " in line.strip():
                g.write(replace_line6+"\n")
                continue

            else:
                g.write(line.strip()+"\n")

def generate_jump_d_param_file(jump_d_params, input_database1, output_prefix, include_contaminants, decoy_generation, decoy_generation_method):
    new_jump_d_params = "jump_d_fly.params"
    fileLines,cometComments,dict1 = storeParamsFile(jump_d_params)
    dict1["input_database1"] = input_database1
    dict1["output_prefix"] = output_prefix
    dict1["include_contaminants"] = include_contaminants
    dict1["decoy_generation"] = decoy_generation
    dict1["decoy_generation_method"] = decoy_generation_method


    for keys in dict1.keys():
        write_log(new_jump_d_params,"{} = {}".format(keys,dict1[keys]))


def generate_QC_param_file(folder, qc_parameter_file, idtxt_list, PSM_recoveray_rate, all_fractions_search):
    try:
        rmFile(qc_parameter_file)
    except:
        print ("  qc_parameter_file is being formed")
    

    run_f_pub = glob.glob(folder+"/*/*.pep.xml")

    if len(run_f_pub) == 0:
        print ("  Comet result is not there so checking for pepXML file")
        run_f_pub = glob.glob(folder+"/*/*.pepXML")

    
    #keep only pepXML files that are in the all_fractions_search list

    run_dict = {}
    for pepfiles in run_f_pub:

        in_string = pepfiles.split(".pep")[0]
        run = pepfiles.split("/")[-2]
        # print (run,"\t", all_fractions_search)
        if run in all_fractions_search:
            run_dict[run] = in_string
#             suffix =  re.search("(\d+)",run).group(1)

#             run_dict[int(suffix)] = in_string

    run_list_sorted = sorted(run_dict.keys())

    all_run_Sorted = []

    for val in run_list_sorted:
        all_run_Sorted.append(run_dict[val])

    for i, files in enumerate(all_run_Sorted):
        write_log (qc_parameter_file,"run"+str(i+1)+" = "+files)


    write_log(qc_parameter_file,"confident_IDtxt = {}".format(idtxt_list[0])) 
    write_log(qc_parameter_file,"accepted_IDtxt = {}".format(idtxt_list[1])) 
    write_log(qc_parameter_file,"PSM_recoveray_rate = {}".format(PSM_recoveray_rate))
    write_log(qc_parameter_file,"output_folder = accepted_PSMs")
                
def prepare_jump_f_input(folder, fdr, pitfile, database, stage, search_engine, min_XCorr):
    
    batch_f_pub = glob.glob(folder+"/*/*.pep.xml")

    if len(batch_f_pub) == 0:
        print ("  Comet result is not there so checking for pepXML file")
        batch_f_pub = glob.glob(folder+"/*/*.pepXML")

    add_string="{}_FDR_{}:".format(stage, fdr)
    for pepfiles in batch_f_pub:
        # out_params = folder+"/jump_fc_stage_{}_FDR_.params".format(stage, fdr)
        out_params = "jump_fc_{}_FDR_.params".format(stage, fdr)

        in_string = pepfiles.split(".pep")[0]
        list_pubs = ""
        add_string+=in_string+"\n"
        list_pubs+=add_string

    fdr_add = "FDR = "+fdr
    pitfile_str = "pit_file = "+pitfile 
    database_str = "database = "+database 
    search_engine_str = "search_engine = "+search_engine
    min_XCorr_str = "min_XCorr = "+min_XCorr

    return out_params, list_pubs, fdr_add, pitfile_str,database_str, search_engine_str, min_XCorr_str




#we need to find stage_0 to get the preprocess results from jump 

def find_stage0_folders_jump(idtxt):
    df = pd.read_csv(idtxt, delimiter=";", skiprows=return_skiprows(idtxt, ";", "Peptide"))

    onerow_df = df.iloc[0:1]
    #;/research_jude/rgs01_jude/groups/penggrp/projects/Alzheimer_BannerSunInstitute/penggrp/proteomics/batch0_pooledSamples/w048/w048.1/w048.16996.1.3.spout
    outfile = onerow_df["Outfile"][0]
    #project path where preprocess files are located
    preprocess = dirname(dirname(dirname(outfile)))
    return preprocess
