import pandas as pd
from os.path import dirname
import re
import os
import argparse
from pyteomics import ms2
from all_JUMPnl_functions import *
from immonium_ion_functions import *
import matplotlib.pyplot as plt
import numpy as np
import pickle
from collections import Counter
import configparser
config = configparser.ConfigParser()


def msg(name=None):
    return '''\n\n\npython mainProgram.py jump_immo_ion_nl_unenrich.params\n\n\n '''

parser = argparse.ArgumentParser(description="Program to find precursor neutral loss or immonium ions", prog="JUMPnl_imm",usage=msg())
parser.add_argument("JUMPnl_imm_param", help="JUMPnl_imm parameter file")

args = parser.parse_args()
curr_path = os.getcwd()

params_file = curr_path+"/"+args.JUMPnl_imm_param


# params_file = sys.argv[1]
# #params_file = "/Users/spoudel1/Desktop/PTM_study/JUMP_validation_Abeta/jump_validator.params"
# #params_file = "../parameterFile/map_comet_jump_fJUMP.params"
config.read(params_file)


ms2File = config["NL"]["ms2File"]
idtxt = config["NL"]["idtxt"]
pepxml = config["NL"]["pepxml"]
stage = config["NL"]["stage"]

# ion_types_str = config["NL"]["ion_types"]

tol = float(config["NL"]["tol"])
output_folder = config["NL"]["output_folder"]
intensity_cut = float(config["NL"]["intensity_cut"])
fixed_noise = float(config["NL"]["fixed_noise"])
precursor_charge_select_on = config["NL"]["precursor_charge_select_on"]
neutral_losses = config["NL"]["neutral_losses"]
immonium_ions = config["NL"]["immonium_ions"]

curr_path = os.getcwd()
newDir = "{}/{}".format(curr_path,output_folder)

mkdir(newDir)

# nl_list = neutral_losses.split(",")
# immonium_ions_list = immonium_ions.split(",")


neutral_losses_dict = {}
immonium_ions_dict = {}

# immonium_ions_dict = {}
# for imm_ion in immonium_ions_list:
# 	ion_mz_list = imm_ion.split(":")
# 	immonium_ions_dict[ion_mz_list[0]] = float(ion_mz_list[1])
if immonium_ions != "0":
	immonium_ions_dict = params_to_dict(immonium_ions)
if neutral_losses != "0":
	neutral_losses_dict = params_to_dict(neutral_losses)


# immonium_ions_dict = {"phosphotyrosine":216.0426, "acetyllysine":143.1179, "acetyllysine_cyc":126.0913} 

print ("\n.... Reading ms2 file {}\n".format(ms2File))
ms2_df = ms2_DF(ms2File)

print (".... Reading complete")
check_variable = os.path.basename(dirname(ms2File))

# pepxml = "/research_jude/rgs01_jude/groups/penggrp/projects/Alzheimer_BannerSunInstitute/penggrp/proteomics/batch0_pooledSamples/PTM_paper_2020/Final_Results/Pho_2020/comprehensive_search/p01/p01.1.pepXML"
jump_modAA_dict, jump_mod_dict, sta_AA = getDynStatModsInfoPepXml(pepxml)

print (".... Reading txt file \n")


pho_data_NR = read_idtxt(idtxt, ms2File, stage,sta_AA, jump_mod_dict)
if ("dma" in neutral_losses_dict.keys()) or ("dmg" in neutral_losses_dict.keys()) or ("mma" in neutral_losses_dict.keys()):
    pho_data_NR = pho_data_NR.loc[pho_data_NR.Peptides.str.contains("R@")]



os.chdir(newDir)


ms2DF_2 = ms2_df.merge(pho_data_NR, how="inner", on = "spectrum")
# if ms2DF_2.shape[0] == 0:
# 	check("check.log", check_variable)

newDF = ms2DF_2.rename(columns={"m/z array":"exp_mz_list","intensity array":"intensity_list"})

print (".. Computing precursor neutral loss\n")

print (".. Using intensity cutoff of {} from the base intensity and remove spectra less than {} intensity".format(intensity_cut, fixed_noise))
cleaning_spectra(newDF, intensity_cut, fixed_noise)
print (".. Spectra cleaned")



if neutral_losses != "0":
	precursor_neutral_loss(newDF, precursor_charge_select_on, neutral_losses_dict, tol, intensity_cut, fixed_noise)
	precursor_nl_count(newDF,"_prec_loss")


if immonium_ions != "0":
	precursor_immonium_loss(newDF, immonium_ions_dict, tol, fixed_noise)

print (".. DONE\n")

# newDF["neutral_loss_results"] = np.where(((newDF.prec_h3po4_loss == 1) | (newDF.prec_hpo3_loss == 1)), "precursor neutral loss present", "no precursor neutral loss")

if newDF.shape[0] != 0:
	# newDF.to_pickle("{}_prec_NL_immIon_output.pkl".format(check_variable))
	newDF.to_csv("{}_prec_NL_immIon_output.txt".format(check_variable), sep="\t")
	# check("check.log", check_variable)


print (".. Making pie chart distribution\n")

if immonium_ions != "0":
    if ("acetyllysine" in immonium_ions_dict.keys()) | ("acetyllysine_cyc" in immonium_ions_dict.keys()):
        if (neutral_losses == "0") & (search_engine.lower() == "comet") & (ptm_type == "acetylation"):
            newDF = newDF.loc[newDF.Peptides.str.contains("K#")]
    precursor_immo_count(newDF,pattern_str="_lowest_tolerance")
    make_pie_plot(newDF,tol,plot_axis="outcome_immonium_ion", title = "All_psms_ImmIon_{}".format(int(tol)))


cpFile("../"+params_file, ".")



if neutral_losses != "0":
    make_pie_plot(newDF,tol,plot_axis="outcome_neutral_loss", title = "All_psms_NL_{}".format(int(tol)))



print (".. COMPLETE\n")

