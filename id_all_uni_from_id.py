import pandas as pd
import numpy as np

def fileToDF(file):
    df = pd.read_csv(file, delimiter="\t")
    return df


def split_outfile(row):
    outfile = row.Outfile
    outfile_split = outfile.split("/")[-1].split(".")
    run = outfile_split[0]
    scan = outfile_split[1]
    charge = outfile_split[3]
    return pd.Series([run,scan,charge])


def make_id_all_txt(psms_DF):
    #first make a column run that all the path
    psms_DF[["Run","Scan#","z"]]=psms_DF["Outfile"].str.split(".",expand=True)
    #just add fraction name in Run#
    psms_DF["Run#"] = psms_DF["Run"].apply(lambda x: x.split("/")[-1])
    #drop Run column
    psms_DF.drop(columns=["Run"], inplace=True)
    
    psms_DF["m/z"] = psms_DF["measuredMH"]



def gen_unique_prot_list(psms_DF):
    psms_DF["group"] = psms_DF.group.astype("int")
    psms_DF["subgroup"] = psms_DF.subgroup.astype("int")
    psms_DF_unique = psms_DF.sort_values(by=["group","subgroup"], ascending = [True, True])[["Peptide","Protein","group","subgroup"]]
    psms_DF_unique = psms_DF_unique.drop_duplicates(subset=["Peptide","group"], keep="first")
    
    #if a peptide maps to 2 different groups
    psms_DF_unique = psms_DF_unique.drop_duplicates(subset=["Peptide"], keep="first")
    
    psms_DF_unique["unique_identifier"] = psms_DF_unique.Peptide +"_"+psms_DF_unique.Protein
    
    unique_pep_protein = list(set(psms_DF_unique.unique_identifier))
    return unique_pep_protein


def count_psms_num(psms_DF):
    psms_DF["XCorr"] = psms_DF.XCorr.astype("float")
    psms_DF_NR = psms_DF.drop_duplicates(subset=["Outfile"], keep="first")
    #count psms per peptide
    psm_count_table = psms_DF_NR.Peptide.value_counts()
    psm_count_DF = pd.DataFrame(psm_count_table).reset_index()
    psm_count_DF.rename(columns={"Peptide":"PSM#","index":"Peptide"}, inplace=True)
    
    return psm_count_DF


def get_bestPSMS(psms_DF):
    psms_DF["XCorr"] = psms_DF.XCorr.astype("float")
    #sort XCorr here
    psms_DF = psms_DF.sort_values(by=["XCorr"], ascending=False)
    psms_DF_NR = psms_DF.drop_duplicates(subset=["Peptide","Protein"], keep="first")
    return psms_DF_NR


def create_new_columns(psms_DF):
    mz_cols = list(psms_DF.columns)
    np_arr = psms_DF.to_numpy()
    run_list = []
    scan_list = []
    charge_list = []
    spectrum_list = []
    key_list = []

    for row in np_arr:
        outfile = row[mz_cols.index("Outfile")]
        ptm_stage = row[mz_cols.index("ptm_stage")]
        outfile_split = outfile.split("/")[-1].split(".")
        run = outfile_split[0]
        scan = outfile_split[1]
        charge = outfile_split[3]
        run_list.append(run)
        scan_list.append(scan)
        charge_list.append(charge)
        spectrum = run+"."+scan+"."+charge
        spectrum_list.append(spectrum)
        key = ptm_stage+"_"+spectrum
        key_list.append(key)
    
    
    psms_DF["Run#"]=run_list
    psms_DF["Scan#"]=scan_list
    psms_DF["z"]=charge_list
    psms_DF["spectrum"]=spectrum_list
    psms_DF["Q-value(%)"] = len(charge_list)*["NA"]
    psms_DF["key"]=key_list
    
 
def generate_id_uni_all_txt_files(pitfile, psms_DF, merge_directory):   
    
	#use merged idtxt to generate id_uni_pep.txt and id_all_pep.txt
	pitDF = fileToDF(pitfile)
	# psms_count df
	psm_no_df = count_psms_num(psms_DF)
	#get best psms
	psms_DF_NR = get_bestPSMS(psms_DF)
	# get unique peptide protein pair
	unique_pep_protein = gen_unique_prot_list(psms_DF)
	#merge psm count list
	psms_DF_NR = psms_DF_NR.merge(psm_no_df, on="Peptide")

	# get required columns for uni_id_pep.txt and uni_all_pep.txt from pit file
	rename_cols_dict = {"SJPGnumber":"Protein Group#","GroupName":"GN","ProteinName":"Protein","FullDescription":"Protein Description"}
	pitDF_extract = pitDF[['SJPGnumber', 'GroupName','ProteinName','FullDescription']]
	pitDF_extract.rename(columns=rename_cols_dict, inplace=True)


	# merge with pit file and non redundant psms matrix
	psms_DF_NR_pitmerged = psms_DF_NR.merge(pitDF_extract, how="inner", on="Protein")
	psms_DF_NR_pitmerged["m/z"] = psms_DF_NR_pitmerged["measuredMH"]

	# create new columns for id_uni/all_pep.txt
	create_new_columns(psms_DF_NR_pitmerged)
	columns_required = ["Peptides","Protein Group#","Protein Accession #","Protein Description","GN","PSM#",
	    "Run#","Scan#","m/z","z", "ppm","XCorr","dCn","Q-value(%)","spectrum", "ptm_stage", "key","Peptides_original",]

	rename_final_cols = {"Peptide":"Peptides","Protein":"Protein Accession #"}
	psms_DF_NR_pitmerged.rename(columns=rename_final_cols, inplace=True)
	id_all_pep_df = psms_DF_NR_pitmerged[columns_required]

	id_all_pep_df["unique_identifier"] = id_all_pep_df.Peptides +"_"+id_all_pep_df["Protein Accession #"]

	id_uni_pep_df = id_all_pep_df.loc[id_all_pep_df["unique_identifier"].isin(unique_pep_protein)]
	id_uni_pep_df_final = id_uni_pep_df[columns_required]

	id_all_pep_df.to_csv(merge_directory+"/publications/id_all_pep.txt", sep="\t", index=None)
	id_uni_pep_df_final.to_csv(merge_directory+"/publications/id_uni_pep.txt", sep="\t", index=None)



