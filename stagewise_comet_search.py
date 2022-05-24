import collections
import sys
import re
import os
import pandas as pd
import numpy as np
from datetime import datetime
import glob

import os.path
import time
import math
from os.path import dirname
import fileinput
from logFunctions import *
from format_pep_xml import *

def return_skiprows(file, delimiter, peptide):
    with open(file, "r") as f:
        skiprows = 0
        for line in f:
            if peptide+delimiter in line:
                break
            else:
                skiprows+=1
    return skiprows


def createOutfile(row, df):
    columns = list(df.columns)
    if "Outfile" in columns:
        outfile_split = row.Outfile.split(".")
        run = outfile_split[-5].split("/")[-1]
        scan = int(outfile_split[-4])
        charge = outfile_split[-2]
        spectrum = run+"."+str(scan)+"."+str(charge)
    else:
        run = row["Run#"]
        scan = row["Scan#"]
        charge = row["z"]
        spectrum = run+"."+str(scan)+"."+str(charge)
    #   print (spectrum)
    return spectrum


#this function converts ms2 to dataframe also used as input to remove accepted psms
def ms2ToDf(ms2File):
    g = open(ms2File,"r")
    lines = g.readlines()
    scan_list = []
    charge_list = []
    MH_list = []
    precursorIon_list = []
    
    ms2_mz = []
    ms2_int = []

    mz_list = []

    for line in lines:
        if "S\t" in line:
            if len(mz_list) >= 1:
                ms2_mz.append(mz_list)
                ms2_int.append(int_list)
            mz_list = []
            int_list = []
    #         print ("spectrum found")
            temp_line = line.strip().split("\t")
            scan = int(temp_line[1])
            scan_list.append(scan)
            precursorIon = float(temp_line[-1])
            precursorIon_list.append(precursorIon)
        elif "Z\t" in line:
    #         print ("Charge found")
            temp_line = line.strip().split("\t")
            charge = int(temp_line[1])
            charge_list.append(charge)

            MH = float(temp_line[-1])
            MH_list.append(MH)


        elif re.match(r"^\d+.*$",line):
    #         print ("ms2 FOund")
            #temp_line = line.strip().split("\t")
            temp_line = re.split("\t| ",line) #to adapt to Zuo-Fei's ms2 pattern
            mz_list.append(float(temp_line[0]))
            int_list.append(float(temp_line[1]))

    ms2_mz.append(mz_list)
    ms2_int.append(int_list)

    dict1 = {"scan":scan_list,"charge":charge_list,"[M+H]+":MH_list,
                "prec_MZ":precursorIon_list,"m/z":ms2_mz,"intensity":ms2_int}

   

    ms2Df = pd.DataFrame.from_dict(dict1)
    
    return ms2Df


#this functons converst dtas to ms2 and it uses the accepted scans dictionary as the tool to removed accepted psms
def dta_to_ms2_select_scans(dtas, new_ms2, scanChargeDf):
    proton = 1.00727646677
    hydrogen = 1.00782503207
    
    f = open(dtas,"r")
    line=f.readlines()
    dtas_dict = {}
    count_outfile = 0
    for x in range(0, len(line),3):
        dta = line[x]
        mass_ms2 = line[x+1]
        ms2_int = line[x+2]
        dta_info = dta.split(".")
        file = dta_info[0]
        scan = dta_info[1]
        
        if scan in list(scanChargeDf.scan):
        
            ppi = dta_info[2]
            charge = dta_info[3]


            neutral_mass_H = float(dta.split()[-2]) #[M+H]
            neutral_mass_only = neutral_mass_H - hydrogen
            neutral_mass_prot = str(neutral_mass_only + proton) #becuase jump has MH so we make it MH+

            check = dta_info[2]
            prec_mz = str(precMZCalc(neutral_mass_only, float(charge)))

            count_outfile+=1
            if int(scan) not in dtas_dict:
                dtas_dict[int(scan)] = [[dta,mass_ms2,ms2_int,dta_info,file,scan, charge,neutral_mass_prot, prec_mz]]
            else:
                dtas_dict[int(scan)].append([dta,mass_ms2,ms2_int,dta_info,file,scan, charge,neutral_mass_prot, prec_mz])
    now = datetime.now()
    print("now =", now)
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%m/%d/%Y %H:%M %p")
    year = now.strftime("%Y")
    date_time = dt_string.split() #date is 0 index, time is 1 index, AM/PM is 2 index

    header_ms2 = "H\tCreationDate\t"+dt_string+"\nH\tExtractor\tMakeMS2\nH\tExtractorVersion\t1.0\nH\tComments\tMakeMS2 written by Suresh Poudel, "+year+"\nH\tExtractorOptions\tMS2/MS1\n"

    #new_ms2 = mzxml_file.split(".mzXML")[0]+".ms2"
    with open(new_ms2,"w") as new_ms2:
        new_ms2.write(header_ms2)
        od = collections.OrderedDict(sorted(dtas_dict.items()))
        count_key = 0
        for key, all_values in od.items():
            count_key+=1
#             if len(all_values) > 1:
#                 print ("\nThis scan has multiple charges, the scan = ", key)
            for mul_values in all_values:
                new_ms2.write("S\t"+str(key)+"\t"+str(key)+"\t"+mul_values[-1]+"\nZ\t"+mul_values[-3]+"\t"+mul_values[-2]+"\n")
                all_mass_list = mul_values[1].split()
                all_int_list = mul_values[2].split()
                for index, val in enumerate(all_mass_list):
                    new_ms2.write(val+"\t"+all_int_list[index]+"\n")
                    
#makes a dataframe with the accepted psms 
#input = idtxt
def idtxt_exclude_acceptedPSMs(idtxt):
    
    with open(idtxt, "r") as f:
        skiprows = 0
        for line in f:
            if "Peptide;" in line:
                break
            else:
                skiprows+=1


    idtxt_df_all = pd.read_csv(idtxt,skiprows=skiprows, delimiter=";")
    idtxt_df = idtxt_df_all.copy()
#     idtxt_df=idtxt_df.loc[idtxt_df.Outfile.str.contains(check_exp)==True]
    idtxt_df["scan"]=idtxt_df.Outfile.apply(lambda x: str(int(x.split(".")[-4])))
    idtxt_df["charge"] = idtxt_df.Outfile.apply(lambda x: x.split(".")[-2])
    scanChargeDf = idtxt_df[["Outfile","scan","charge"]]
    
    return scanChargeDf


def idtxt_exclude_acceptedPSMs_acceptedProtein(idtxt):
    
    with open(idtxt, "r") as f:
        skiprows = 0
        for line in f:
            if "Peptide;" in line:
                break
            else:
                skiprows+=1


    idtxt_df_all = pd.read_csv(idtxt,skiprows=skiprows, delimiter=";")
    idtxt_df = idtxt_df_all.copy()
#     idtxt_df=idtxt_df.loc[idtxt_df.Outfile.str.contains(check_exp)==True]
    idtxt_df["scan"]=idtxt_df.Outfile.apply(lambda x: str(int(x.split(".")[-4])))
    idtxt_df["charge"] = idtxt_df.Outfile.apply(lambda x: x.split(".")[-2])
    scanChargeDf = idtxt_df[["Outfile","scan","charge"]]
    
    accepted_protein_list = list(set(idtxt_df["Protein"]))

    return accepted_protein_list,scanChargeDf


def additional_qc_parse(additional_qc_file, logFile):
    df = pd.read_csv(additional_qc_file, header=None)

    np_arr = df.to_numpy()
    idtxt = np_arr[0][0]
    accepted_protein_list,scanChargeDf=idtxt_exclude_acceptedPSMs_acceptedProtein(idtxt)
    write_log (logFile,"  Total scans considered for exclusion from {} = {}".format(idtxt, scanChargeDf.shape[0]))
    for row in np_arr[1:]:
        idtxt2 = row[0]
        accepted_protein_list2,scanChargeDf2=idtxt_exclude_acceptedPSMs_acceptedProtein(idtxt2)
        write_log (logFile,"  Total scans considered for exclusion from {} = {}".format(idtxt2, scanChargeDf2.shape[0]))

        scanChargeDf.append(scanChargeDf2)
        frames = [scanChargeDf, scanChargeDf2]
        scanChargeDf = pd.concat(frames)
    return scanChargeDf

#makes the ms2 from the dtaframe 
def createMS2_fromDF(df,new_folder, check_exp, logFile):
    
    mz_cols = list(df.columns)
    proton = 1.00727646677
    write_log (logFile,"  Generating .ms2 files\n")
    now = datetime.now()
    #print (prev_mod)
    write_log(logFile,"  now =", now)
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%m/%d/%Y %H:%M %p")
    year = now.strftime("%Y")
    date_time = dt_string.split() #date is 0 index, time is 1 index, AM/PM is 2 index
    
    header_ms2 = "H\tCreationDate\t"+dt_string+"\nH\tExtractor\tMakeMS2\nH\tExtractorVersion\t1.0\nH\tComments\tMakeMS2 written by Suresh Poudel, "+year+"\nH\tExtractorOptions\tMS2/MS1\n"
    new_ms2_file = new_folder+"/"+check_exp+".ms2"
    # new_ms2_file = specLibFolder+"/Spectral_Library_Decoy.spLib" #modified the extension name to ms2pep
    write_log (logFile,"  ms2 file is being created")
    
    with open(new_ms2_file,"w") as new_ms2:
        new_ms2.write(header_ms2)

        for row in df.itertuples():
            mzIntDict = {}
            scan = str(row[mz_cols.index("scan")+1])
            massNeutral = float(row[mz_cols.index("[M+H]+")+1])

            precursorMZ = float(row[mz_cols.index("prec_MZ")+1]) #because it is decoy precursor ions using precursor swap
            

            mz = row[mz_cols.index("m/z")+1]
            intensity = row[mz_cols.index("intensity")+1]


            charge = int(row[mz_cols.index("charge")+1])
 


            new_ms2.write("S\t"+scan+"\t"+scan+"\t"+str(precursorMZ)+"\n")
            new_ms2.write("Z\t"+str(charge)+"\t"+str(massNeutral)+"\n")

            for index, val in enumerate(mz):
                new_ms2.write(str(val)+"\t"+str(intensity[index])+"\n")

    write_log (logFile,"  Done creating {}\n".format(new_ms2))


#fucntion to remove accepted psms and generate dataframe to creat ms2 file using the scan that were not accepted
def make_ms2_removed_accepted_psms(scanChargeDf, ms2, new_folder, logFile):

    check_exp = ms2.split("/")[-2]
    scanChargeDf["key"] = scanChargeDf.scan+"\t"+scanChargeDf.charge
    
    scanChargeDf=scanChargeDf.loc[scanChargeDf.Outfile.str.contains(check_exp)==True]
    exclusion_scan_charge =list(set(list(scanChargeDf["key"])))
    
    write_log (logFile,"    QC leads to a total of {} scans removal for fraction {}".format(len(exclusion_scan_charge),check_exp))
    
    #makes dataframe of the entire ms2 file
    ms2_df = ms2ToDf(ms2)
    
    ms2_df["key"] = ms2_df.scan.astype("str")+"\t"+ms2_df.charge.astype("str")
    
    #excludes accepted psms
    ms2_df = ms2_df.loc[~ms2_df.key.isin(exclusion_scan_charge)]
    
    #creates the ms2 from dataframe ms2_df after removal of accepted psms
    createMS2_fromDF(ms2_df,new_folder, check_exp, logFile)
    

    
############## Main function for stagewise searching ###############


#In this function the input files are the mzxmls for each stages
#the function generates fractionwise folders inside stage folders
#it softlinks mzxml files inside those fraction folders
#It also makes preprocessing folder in the Pipeline results folder
#Perform preprocessing for each mzXML files and convert them to ms2 files (mass corrected and preprocessed)
#If the stage = Stage_0 it performs whole proteome searches
#it uses scan\tcharge containing dataframe that have excluded the scans and is used to generate the reduced ms2 files


#folders = stage wise folders
#mzXMLs = list of mzXMLfiles
#scanChargeDf = dataframe generated from function idtxt_exclude_acceptedPSMs(idtxt) uses idtxt as input


def rename_pep_xml(mzXMLs, folders):
    for mz_file in mzXMLs:
        
        basefile = os.path.basename(mz_file)
        basefile_noext = basefile.split(".")[0]
        
        os.chdir(folders+"/"+basefile_noext)
        
        if exists(basefile_noext+".pep.xml"):
            searchPep = basefile_noext+".pep.xml"
        else:
            searchPep = glob.glob(basefile_noext+"*-*.pep.xml")[0]

        makedirectory("comet")
        addSuffixPepMvFiles(searchPep)


def create_job_file(mzxml_file, comet_params, folders, basefile_noext, comet):
    log_dir = folders+"/"+basefile_noext+"/log"
    cmd1 = "rm -r "+log_dir
    try:
        os.system(cmd1)
    except:
        print ("No log directory for ", basefile_noext)
    makedirectory(log_dir)
    
    job_header = "#!/bin/bash\n#BSUB -P TestComet\n#BSUB -J comet\n#BSUB -oo "+log_dir+"/log.out\n#BSUB -eo "+log_dir+"/error.err\n#BSUB -n 8\n"

    job_body1 = comet+" -P"+comet_params+" "+mzxml_file

    #job_body1 = comet+" -P"+comet_params+" "+mzxml_file

    jobfile = basefile_noext+".sh"
    with open(jobfile,"w") as new_file:
        new_file.write(job_header+job_body1)
    return jobfile



def submit_job(jobf,queue,mem):
    cmd = 'bsub -q '+queue+' -R "rusage[mem='+mem+']" < '+jobf
    os.system(cmd)
    

#cnt = 0
#jobs = 4
#folder = where you have subfolders for each fraction after search
#for comet search cnt = 0 jobs = 1
#for preprocess cnt = 0 jobs = 4
def wait_Function(jobs, folder, mzXMLs):
    cnt=0
    while len(mzXMLs)*jobs != cnt:
        cnt = 0
        result = False

        logfiles = glob.glob(folder+"/*/log/log*.out")

        if len(logfiles) != 0:
            for logfile in logfiles:
                result = read_log_file(logfile)
                if result == True:
                    cnt+=1

        time.sleep(30)
        
def select_max_pep_xml(pepxml_list):
    max_suffix = 0
    for index,peps in enumerate(pepxml_list):
        num_suffix_withbase = peps.split(".pep.xml")[0]
        num_suffix = int(num_suffix_withbase.split(".")[-1])
        # print ("\n\n{}\n\n".format(num_suffix))
        if num_suffix > max_suffix:
            max_suffix = num_suffix
            max_suffix_pep = peps 
    # print (max_suffix_pep)
    return max_suffix_pep


def run_tag_program(mzXMLs, folders, jump_tag_program, tags_input_path, cluster):
    for mz_file in mzXMLs:
        os.chdir(folders)
        basefile = os.path.basename(mz_file)
        
        basefile_noext = basefile.split(".")[0]
        tag_file = glob.glob(tags_input_path+"/"+basefile_noext+"*.tags")[0]

        os.chdir(folders+"/"+basefile_noext)
        #this is renamed pepxml file
        searchPep_list = glob.glob(basefile_noext+"*.pep.xml")
        searchPep = select_max_pep_xml(searchPep_list)
        #input for tags is ready 
        #1. searchPep 2. tag_file

        #prepare jump_tag_params file
        generate_tag_match_param_file("tag_qc.params", folders+"/"+basefile_noext+"/"+searchPep, tag_file)

        #submit the job to the server now
        # bsub  "jump -q jump_qc_HH_tmt16_human.params"
        if cluster == "0":
            cmd = "{} {}".format(jump_tag_program,"tag_qc.params")
           
            os.system(cmd)
        else:
            #submit the job to the server now
            # bsub  "jump -q jump_qc_HH_tmt16_human.params"
            cmd = 'bsub -P qc_program -R "rusage[mem=10G]" "{} {}"'.format(jump_tag_program,"tag_qc.params")
            os.system(cmd)

    if cluster == "1":
        wait_function_after_tag_match(mzXMLs, folders, "tag_qc.params")

def run_tag_program_server(mzXMLs, folders, jump_tag_program, dtasFolder):
    for mz_file in mzXMLs:
        os.chdir(folders)
        basefile = os.path.basename(mz_file)
        
        basefile_noext = basefile.split(".")[0]
        tag_file = glob.glob(dtasFolder+"/"+basefile_noext+"*.tags")[0]
        

        os.chdir(folders+"/"+basefile_noext)
        #this is renamed pepxml file
        searchPep = glob.glob(basefile_noext+"*.pep.xml")[0]

        #input for tags is ready 
        #1. searchPep 2. tag_file

        #prepare jump_tag_params file
        generate_tag_match_param_file("tag_qc.params", folders+"/"+basefile_noext+"/"+searchPep, tag_file)

        #submit the job to the server now
        # bsub  "jump -q jump_qc_HH_tmt16_human.params"
        cmd = 'bsub -P qc_program -R "rusage[mem=10G]" "{} {}"'.format(jump_tag_program,"tag_qc.params")
        print (cmd)
        os.system(cmd)


def wait_function_after_tag_match(mzXMLs, folders, tag_qc_params):
    cnt=0
    tag_files = []
    # print (mzXMLs)

    while len(tag_files)!= len(mzXMLs):
        cnt = 0
        result = False

        tag_files = glob.glob(folders+"/*/*/"+tag_qc_params)
        # print (tag_files)
        time.sleep(30)


def stage_wise_search(folders, mzXMLs, comet,logFile,cluster,scanChargeDf=None):
    
    stage = os.path.basename(folders)  
    write_log (logFile,"\n\nWorking on the folder structure for {}".format(stage))
    
    #input mzXMLs file list
    for mz_file in mzXMLs:
        os.chdir(folders)
        comet_params = glob.glob(folders+"/stage_*_comet.params")[0]
        
        basefile = os.path.basename(mz_file)
        basefile_noext = basefile.split(".")[0]
        makedirectory(basefile_noext)
        
        mzxml_file = mz_file.split(".ms2")[0]+".mzXML"

        #copy mzxml file to the stage folder
        softlink_mzxml([mz_file], basefile_noext)

        os.chdir(folders+"/"+basefile_noext)
        cpFile(mz_file, folders+"/"+basefile_noext)
        cpFile(comet_params, "comet.params")
        
        if cluster == "0":
            #comet search command
            job_body1 = "{} -P{} {}".format(comet, "comet.params", "*.ms2 > search_log.txt")
            #run comet
            os.system(job_body1)
        else:
            jobfile = create_job_file("*.ms2", "comet.params", folders, basefile_noext, comet)
            submit_job(jobfile,"standard","20G")
            write_log (logFile,"  Comet search for {} is submitted".format(mz_file)) 
    
    if cluster == "1":
        wait_Function(1, folders, mzXMLs)

def dta_to_ms2(dtas, new_ms2):
    proton = 1.00727646677
    hydrogen = 1.00782503207
    f = open(dtas,"r")
    line=f.readlines()
    dtas_dict = {}
    count_outfile = 0
    for x in range(0, len(line),3):
        dta = line[x]
        mass_ms2 = line[x+1]
        ms2_int = line[x+2]
        dta_info = dta.split(".")
        file = dta_info[0]
        scan = dta_info[1]
        ppi = dta_info[2]
        charge = dta_info[3]
       

        neutral_mass_H = float(dta.split()[-2]) #[M+H]
        neutral_mass_only = neutral_mass_H - hydrogen
        neutral_mass_prot = str(neutral_mass_only + proton) #becuase jump has MH so we make it MH+

        check = dta_info[2]
        prec_mz = str(precMZCalc(neutral_mass_only, float(charge)))

        count_outfile+=1
        check_key = scan+"\t"+charge
        
        if int(scan) not in dtas_dict:
            dtas_dict[int(scan)] = [[dta,mass_ms2,ms2_int,dta_info,file,scan, charge,neutral_mass_prot, prec_mz]]
        else:
            dtas_dict[int(scan)].append([dta,mass_ms2,ms2_int,dta_info,file,scan, charge,neutral_mass_prot, prec_mz])
    now = datetime.now()
    print("now =", now)
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%m/%d/%Y %H:%M %p")
    year = now.strftime("%Y")
    date_time = dt_string.split() #date is 0 index, time is 1 index, AM/PM is 2 index

    header_ms2 = "H\tCreationDate\t"+dt_string+"\nH\tExtractor\tMakeMS2\nH\tExtractorVersion\t1.0\nH\tComments\tMakeMS2 written by Suresh Poudel, "+year+"\nH\tExtractorOptions\tMS2/MS1\n"

    #new_ms2 = mzxml_file.split(".mzXML")[0]+".ms2"
    with open(new_ms2,"w") as new_ms2:
        new_ms2.write(header_ms2)
        od = collections.OrderedDict(sorted(dtas_dict.items()))
        count_key = 0
        for key, all_values in od.items():
            count_key+=1
#             if len(all_values) > 1:
#                 print ("\nThis scan has multiple charges, the scan = ", key)
            for mul_values in all_values:
                new_ms2.write("S\t"+str(key)+"\t"+str(key)+"\t"+mul_values[-1]+"\nZ\t"+mul_values[-3]+"\t"+mul_values[-2]+"\n")
                all_mass_list = mul_values[1].split()
                all_int_list = mul_values[2].split()
                for index, val in enumerate(all_mass_list):
                    new_ms2.write(val+"\t"+all_int_list[index]+"\n")
              
    #print ("\nTotal dta keys = ",count_key)
#jump -deisotope jump_ss_HH_tmt10_mouse.params mix_ratio.mzXML


# def precMZCalc(MH, z): #MH = MH+ from dta file, z = charge and proton = H+
#     proton = 1.00727646677
#     hydrogen = 1.00782503207
#     precmz = (float(MH)+((int(z)-1)*proton))/int(z)
#     return precmz


def precMZCalc(M, z): #M = Neutral mass
    proton = 1.00727646677
    precmz = (M+(z*proton))/z
    return precmz


def generate_tag_match_param_file(tag_qc_params, pepxml, tagfile):
    try:
        rmFile(tag_qc_params)
    except:
        print ("  tag_qc_params is being formed")
    
    write_log(tag_qc_params, "[tagqc]")
    write_log(tag_qc_params,"pep_xml = {}".format(pepxml)) 
    write_log(tag_qc_params,"tag_file = {}".format(tagfile)) 
    write_log(tag_qc_params,"start_scan = 0")
    write_log(tag_qc_params,"end_scan = 0")
    write_log(tag_qc_params,"min_tag_length = 2")
