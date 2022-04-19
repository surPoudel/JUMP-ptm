import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pepxml_generation_reorder import write_log


def make_scores_spectrumDF(dfMz_valid):
    weight_Eval = []
    peptide_list = []
    xcorr_list = []
    expect_list = []
    expect_score_list = []
    protein_list = []
    spectrum_list = []
    tag_score_list = []
    tag_list = []
    rank_list = []
    mod_peptide_list = []
#     expect_zscore_list=[]
#     tag_zscore_list = []

    mz_cols = list(dfMz_valid.columns)
    np_arr = dfMz_valid.to_numpy()
    for row in np_arr:
        spectrum = str(row[mz_cols.index("spectrum_jump")])
        search_hits = row[mz_cols.index("search_hit_reordered")]


        for index,val in enumerate(search_hits):
            rank_list.append(val["hit_rank"])
            jumpTag = ""
            jumpTagScore = 0
            weight_Eval.append(val["search_score"]['WeightedEvalue'])
            peptide_list.append(val['peptide'])
            mod_peptide_list.append(val["modified_peptide"])
            xcorr_list.append(val["search_score"]['xcorr'])
            expect = val["search_score"]['expect']
            expect_list.append(expect)
            expect_score_list.append(-1*(np.log10(expect)))
#             expect_z_score = zscore(val["search_score"]['expect'], comet_eval_mean, comet_eval_std)
#             expect_zscore_list.append(expect_z_score)
            try:
                jumpTagScore = val["search_score"]['jumpTagScore']
                jumpTag = val["search_score"]['jumpTag']
            except:
                pass
            tag_score_list.append(jumpTagScore)
            tag_list.append(jumpTag)
#             tag_z_score = zscore(tag_score_list, tag_mean, tag_std)
#             tag_zscore_list.append(tag_z_score)
            protein_list.append(val['proteins'][0]['protein'])
            spectrum_list.append(spectrum)
#     df_scores = pd.DataFrame([spectrum_list,protein_list, peptide_list,weight_Eval,xcorr_list,expect_list,tag_list,tag_score_list,expect_zscore_list,tag_zscore_list]).T
#     df_scores.columns = ["spectrum","proteins","peptide","weightedEvalue","xcorr","expect","jump_Tag","jumpTagScore","expect_zscore","tag_zscore"]
    
    df_scores = pd.DataFrame([spectrum_list,rank_list,protein_list, peptide_list,mod_peptide_list,weight_Eval,xcorr_list,expect_list,expect_score_list,tag_list,tag_score_list]).T
    df_scores.columns = ["spectrum","rank","proteins","peptide","mod_peptide","weightedEvalue","xcorr","expect","-log10(expect)","jump_Tag","jumpTagScore"]
    df_scores["Type"] = df_scores.apply(typePeptide, axis=1)
    return df_scores



###Target decoy test ###
def countcumulativeTarget(row):
    if row["Type"] == "Target":  #Define Target 
        value = 1
    else:
        value =0
    return value

def countcumulativeDecoy(row):
    if row["Type"] == "Decoy":  #Define Decoy
        value = 1
    else:
        value =0
    return value



def calcFDR(row):
    if row.cumsumTarget == 0:
        FDR = 100
    else:
        FDR = row.cumsumDecoy/row.cumsumTarget*100
    return FDR



def FDR_Target_Decoy(df,logFile,sortCol="JDscore"):
    reqd_cols = list(df.columns)
    df1 = df.copy()
    #define Target Decoy
    # df1["Target-Decoy"]=df1.Peptide_ID.apply(lambda row: "Decoy" if "Decoy" in row else "Target")
    if sortCol != "expect":
        df2=df1.sort_values([sortCol,"Type"],ascending=[False,False])   #for labeled dataset
    else:
        df2=df1.sort_values([sortCol,"Type"],ascending=[True,False])   #for labeled dataset
    df2["Target_value"] = df2.apply(countcumulativeTarget, axis=1)
    df2["Decoy_value"] = df2.apply(countcumulativeDecoy, axis=1)
    #df['SUM_C'].cumsum()
    df2["cumsumTarget"] = df2["Target_value"].cumsum()
    df2["cumsumDecoy"] = df2["Decoy_value"].cumsum()
    df2["FDR"] = df2.apply(calcFDR, axis=1)
    if "FDR" in reqd_cols:
        addedColsAll = reqd_cols
    else:
        addedColsAll = reqd_cols+["FDR"]

    # print ("    Total Targets and Decoys for {} is {}".format(sortCol,df2.Type.value_counts()))
    write_log (logFile,"    Total Targets and Decoys for {} is {}".format(sortCol,df2.Type.value_counts()))
    return df2[addedColsAll]


def histogramPlot(matched_df, unmatched_df, xaxis, figname,label1, label2): #bins2 for xaxis label
    minv = np.min(unmatched_df[xaxis])
    maxv = np.max(matched_df[xaxis])
    bins = np.linspace(minv,maxv)

    plt.rcParams.update({'font.size': 10})
    fig,ax = plt.subplots(figsize=(4,2))
    plt.yticks(color="black")
    # size, scale = 1000, 10

    commutes2 = matched_df[xaxis]
    commutes2.plot.hist(grid=False, bins=bins, rwidth=0.9,
                   color='#F4F6F7',edgecolor='black', linewidth=1.0)
    commutes = unmatched_df[xaxis]
    commutes.plot.hist(grid=False, bins=bins, rwidth=0.9,
                   color='#808B96',edgecolor='black', linewidth=1.0)
    # bins_labels(bins2)
    # Hide grid lines
    # plt.grid(False)

    plt.title('')
    plt.xlabel(xaxis)
    plt.ylabel('Number of PSMs')
#     plt.xlim(xlim)
    #   plt.xticks(bins)
    # plt.grid(axis='y', alpha=0.75)
    plt.legend([label1, label2],loc="best")
    figurename = figname+".pdf"
    figurename1 = figname+".png"
    fig.savefig(figurename, bbox_inches="tight", dpi=600 )
    fig.savefig(figurename1, bbox_inches="tight", dpi=600 )

def typePeptide(row):
    peptide = row.proteins
    if "Decoy" in peptide:
        pepType = "Decoy"
    else:
        pepType = "Target"
        
    return pepType

####
#combined evalue with tag information only
####orignial tag only based separation
def qc_target_decoy(df, logFile,search_hits_column = "search_hit_reordered", separation="all", targe_decoy_score = "WeightedEvalue"):
    searchhits_columns_xcorr = list(df[search_hits_column])
    spectrum_list_all = list(df["spectrum"])
    xcorr_list = []
    protein_list_xcorr = []
    spectrum_list = []
    
    for index, x in enumerate(searchhits_columns_xcorr):
#         print (index)
        for val in x:
            if val["hit_rank"] == 1:
                
                if separation == "tagOnly":
                    
                    if 'jumpTag' not in list(val['search_score'].keys()):
                        pass
                    else:
                        protein_list_xcorr.append(val['proteins'][0]['protein'])
                        xcorr_list.append(val['search_score'][targe_decoy_score])
                        spectrum_list.append(spectrum_list_all[index])
                        
                else:
                    protein_list_xcorr.append(val['proteins'][0]['protein'])
                    xcorr_list.append(val['search_score'][targe_decoy_score])
                    spectrum_list.append(spectrum_list_all[index])

        
    df_xcorr = pd.DataFrame([spectrum_list,protein_list_xcorr,xcorr_list]).T
    df_xcorr.columns = ["spectrum","proteins",targe_decoy_score]
    if targe_decoy_score == "expect":
        df_xcorr["-log10(eval)"] = df_xcorr.expect.apply(lambda x: -1*(np.log10(x)))
        
        targe_decoy_score = "-log10(eval)"
    df_xcorr["Type"] = df_xcorr.apply(typePeptide, axis=1)
#     print (df_xcorr)
    df2=FDR_Target_Decoy(df_xcorr,logFile,targe_decoy_score)
    df3 = df2.query('FDR<1')
    
    write_log (logFile,"Total ID at FDR = 1% at psms level {}".format(df3.shape[0]))
    # print ("Total ID at FDR = 1% at psms level {}".format(df3.shape[0]))
    
    return df2
    




def qc_figures(df_xcorr, scoreCol, figname):
    target = df_xcorr.loc[df_xcorr.Type == "Target"]
    decoy = df_xcorr.loc[df_xcorr.Type == "Decoy"]
    histogramPlot(target, decoy, scoreCol, figname,"target", "decoy") #bins2 for xaxis label
    


def histogramPlot_zscores(matched_df, xaxis, figname,label1): #bins2 for xaxis label
    minv = np.min(matched_df[xaxis])
    maxv = np.max(matched_df[xaxis])
    bins = np.linspace(minv,maxv)

    plt.rcParams.update({'font.size': 10})
    fig,ax = plt.subplots(figsize=(4,2))
    plt.yticks(color="black")
    # size, scale = 1000, 10

    commutes2 = matched_df[xaxis]
    commutes2.plot.hist(grid=False, bins=bins, rwidth=0.9,
                   color='#F4F6F7',edgecolor='black', linewidth=1.0)


    plt.title('')
    plt.xlabel(xaxis)
    plt.ylabel('Number of PSMs')
#     plt.xlim(xlim)
    #   plt.xticks(bins)
    # plt.grid(axis='y', alpha=0.75)
    plt.legend([label1],loc="best")
    figurename = figname+".pdf"
    figurename1 = figname+".png"
    fig.savefig(figurename, bbox_inches="tight", dpi=600 )
    fig.savefig(figurename1, bbox_inches="tight", dpi=600 )

