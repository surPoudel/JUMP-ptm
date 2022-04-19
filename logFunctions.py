import pandas as pd
import os, sys
import glob

def write_log(logFile,*args):
    with open(logFile, "a") as log_file:
        line = ' '.join([str(a) for a in args])
        log_file.write(line+'\n')
        print(line)
        
def rmFile(filename):
    cmd1 = "rm "+filename
    os.system(cmd1)

def cpFile(filename, destinationFile, mode="cp"):
    cmd1 = "{} {} {}".format(mode, filename, destinationFile)
    os.system(cmd1)



def makedirectory(folName):
    #create search output directory
    cmdDir = "mkdir "+folName
    
    try:
        os.system(cmdDir)
    except:
        print ("Directory exist")

def fileToDF(file):
    df = pd.read_csv(file, delimiter="\t")
    return df


def exists(path):
    """Test whether a path exists.  Returns False for broken symbolic links"""
    try:
        st = os.stat(path)
    except os.error:
        return False
    return "The {} file exists. The program will continue".format(path)

def softlink_mzxml(list_of_mzxml, destination_folder):
    for mzfile in list_of_mzxml:
        #ln -s /home/spoudel1/PanPTM_Paper_2021/Pipeline/*.mzXML Stage_1
        cmd1 = "ln -s {} {}".format(mzfile, destination_folder)
        os.system(cmd1)
        
        
def addSuffixPepMvFiles(pepxml):
    total_pep_xml = glob.glob("comet/*.pep.xml")
    
    suffix = "{}".format(len(total_pep_xml)+1)
    filesplit = pepxml.split(".")
    txt_file = pepxml.split(".pep.xml")[0]+".txt"
    
    new_pep_xml = filesplit[0]+".{}.pep.xml".format(suffix)
    new_txt_File = filesplit[0]+".{}.txt".format(suffix)
    
    cpFile(pepxml, new_pep_xml, mode="cp")
    cpFile(txt_file, new_txt_File, mode="cp")
    
    cpFile(pepxml, "comet", mode="mv")
    cpFile(txt_file, "comet", mode="mv")
    
    
def read_log_file(logfile):
    result = False
    with open(logfile, 'r') as f:
        all_lines = f.readlines()
        if len(all_lines) > 2:
            second_last_line = all_lines[-2]
            last_line = all_lines[-1]

            if ("for stderr output of this job" in last_line) or ("for stderr output of this job" in second_last_line):
                result = True
        else:
            result = False
    return result


