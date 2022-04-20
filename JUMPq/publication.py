import os
import pandas as pd

def return_skiprows(file, delimiter, peptide):
    with open(file, "r") as f:
        skiprows = 0
        for line in f:
            if peptide+delimiter in line:
                break
            else:
                skiprows+=1
    return skiprows

def makedirectory(folName):
    #create search output directory
    cmdDir = "mkdir "+folName
    
    try:
        os.system(cmdDir)
    except:
        print ("Directory exist")


def generateTables(dfPep, dfProt, params):
    dirPub = os.path.dirname(params['idtxt']) + "/publications/"
    print (dirPub)
    ########################
    # Peptide-level tables #
    ########################
    dfPep = dfPep.reset_index()
    dfPep = dfPep.rename({"index": "Peptides"}, axis=1)

    # Unique peptides
    unique_id_pep_file = os.path.join(dirPub, "id_uni_pep.txt")
    unique_all_pep_file = os.path.join(dirPub, "id_all_pep.txt")
    dfJumpf = pd.read_table(unique_id_pep_file, sep="\t", skiprows=return_skiprows(unique_id_pep_file, "\t", "Peptides"), header=0)

    # print ("dfJumpf columns",dfJumpf.columns)
    # print ("dfPep columns",dfPep.columns)

    dfUniPep = dfJumpf.merge(dfPep, on="Peptides")

    # All peptides
    dfJumpf = pd.read_table(unique_all_pep_file, sep="\t", skiprows=return_skiprows(unique_all_pep_file, "\t", "Peptides"), header=0)
    dfAllPep = dfJumpf.merge(dfPep, on="Peptides")

    ########################
    # Protein-level tables #
    ########################
    # dfProt = dfProt.reset_index()
    # dfProt = dfProt.rename({"index": "Protein Accession #"}, axis=1)

    # unique_id_prot_file = os.path.join(dirPub, "id_uni_prot.txt")
    # unique_all_prot_file = os.path.join(dirPub, "id_all_prot.txt")

    # # Unique proteins
    # dfJumpf = pd.read_table(unique_id_prot_file, sep="\t", skiprows=return_skiprows(unique_id_prot_file, "\t", "Peptides"), header=0)
    # dfUniProt = dfJumpf.merge(dfProt, on="Protein Accession #")

    # # All peptides
    # dfJumpf = pd.read_table(unique_all_prot_file, sep="\t", skiprows=return_skiprows(unique_all_prot_file, "\t", "Peptides"), header=0)
    # dfAllProt = dfJumpf.merge(dfProt, on="Protein Accession #")

    # return dfUniPep, dfAllPep, dfUniProt, dfAllProt
    return dfUniPep, dfAllPep