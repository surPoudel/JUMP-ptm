
import pandas as pd
from pyteomics import mass
import numpy as np

#this function is different in consensus library as we cannot sum modifications there
#as the PTMs are used to extract unimod information
def spectrumToDict(spectrum):
    dict1 = {}
    spectrumCommaSplit = spectrum.split(",")
    for x in spectrumCommaSplit:
        y=x.split("_")
        if y[0] not in dict1.keys():
            dict1[y[0]] = [float(y[2])]
        else:
            dict1[y[0]].append(float(y[2]))
    dict2 = {}

    for key in dict1.keys():
        value = np.sum(list(dict1[key]))
        # value = np.sum(list(set(dict1[key])))
        dict2[key] = value
    return dict2

def ionSeriesIonLossSpeRes(peptide,massPosDict, maxcharge=1,useMod ="Yes"):
    #first check if the mods has phosphorylation if yes get the maximum and minimum position
    phoPos = getPhoPosition(massPosDict)
    minPos = 1    #assigns minPosition which is later updated
    maxPos = len(peptide) #assigns maximum position which is later updated
    if len(phoPos) >= 1:
        minPos = np.min(phoPos) #computes minimum phospho position so tha b and y ions can be computed according for phospho loss
        maxPos = np.max(phoPos) #computes maximum phospho position so that b and y ions can be computed according ly

    if useMod != "Yes": #this checks whether the modificaiton is searched or not, generally this is always Yes for localization 
        massPosDict = {0:0.0} #if modificaiton is no than the dictionary has no information
    h2o = mass.calculate_mass(formula='H2O') #h2o mono mass
    co = mass.calculate_mass(formula='CO') #co mono mass
    nh3 = mass.calculate_mass(formula='NH3') #nh3 mono mass
    xmassMore = co-mass.calculate_mass(formula='H2') #hydrogen mono mass
    proton = mass.calculate_mass(formula='H+') #proron mono mass
    hydrogenAtom = mass.calculate_mass(formula='H') #hydrogen atom mono mass
    hpo3 = mass.calculate_mass(formula='HPO3') #computes hpo3 mono mass
    h3po4 = mass.calculate_mass(formula='H3PO4') #computes h3po4 mono mass

    all_ions_dict = {}#iniitates a dictionary for theoretical ions which is later conveted to dataframe
#possible ion loss or no ion loss
    ionloss = {"":0,"-H2O":h2o,"-HPO3":hpo3,"-H3PO4":h3po4,"-NH3":nh3}
    
#this section computes a,b,c ions    
    for i in range(1, len(peptide)+1):
        addModMass = 0.0
        for massKey in massPosDict.keys():
            if int(massKey) <= i:
                addModMass += float(massPosDict[massKey])
                
        valAddKey(all_ions_dict,"Seq",peptide[i-1])
#         print (peptide[0:i-1])
        for losses in ionloss.keys():
            if losses == "-H2O":
#water loss if the amino acid sequence is STED
                fate = checkAA(peptide[0:i],["S","T","E","D"])
                for charge in range(1, maxcharge):
                    if (fate == "True") and (i < len(peptide)):
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
                    else:
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
            if losses == "-NH3":
#ammonia loss if the amono acid sequence is RKQN
                fate = checkAA(peptide[0:i],["R","K","Q","N"])
                for charge in range(1, maxcharge):
                    if (fate == "True") and (i < len(peptide)):
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
                        valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
                    else:
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),"---")
                                            
            if (losses == "-H3PO4") or (losses == "-HPO3"):
                for charge in range(1, maxcharge):
                    if "79.96" in str(massPosDict.values()):
#                         fate = checkAA(peptide[0:i],["S","T","Y"])
#                         if fate == "True":
#phosoho loss if massPosDict have phopsho loss in it. If yes based on minPos and maxPos ions are calculated
                        if (i >= minPos) and (i < len(peptide)):
                            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
                            valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
                        else:
                            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
                            valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),"---")
                                            
            if losses == "":
                for charge in range(1, maxcharge):
                    if i < len(peptide):
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
                        valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
                        valAddKey(all_ions_dict,"c"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge)
                        valAddKey(all_ions_dict,"c(-1)"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge-hydrogenAtom)
                        valAddKey(all_ions_dict,"c1"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge+hydrogenAtom)
                        valAddKey(all_ions_dict,"c2"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge+(2*hydrogenAtom))

                    else:
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"c"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"c(-1)"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"c1"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"c2"+losses+"+"+str(charge),"---")

#this section computes x,y,z ions
    for i in range(0, len(peptide)):
        
        addModMass = 0.0
        for massKey in massPosDict.keys():
            if int(massKey) > i:
                addModMass += float(massPosDict[massKey])
        for losses in ionloss.keys():    
            if losses == "-H2O":
                fate = checkAA(peptide[i:len(peptide)],["S","T","E","D"])
                for charge in range(1, maxcharge):
                    if (fate == "True") and (i >=1):
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
                    else:
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
            
            if losses == "-NH3":
                fate = checkAA(peptide[i:len(peptide)],["R","K","Q","N"])
                for charge in range(1, maxcharge):
                    if (fate == "True") and (i >=1):
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
                    else:
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
            
            
            if (losses == "-H3PO4") or (losses == "-HPO3"):
#                 fate = checkAA(peptide[i:len(peptide)],["S","T","Y"])
                for charge in range(1, maxcharge):
                    if "79.96" in str(massPosDict.values()):
                        if (i < maxPos) and (i >=1):
                            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
                        else:
                            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
            if losses == "":
                for charge in range(1, maxcharge):
                    if (i >=1):
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
                        valAddKey(all_ions_dict,"x"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)+xmassMore-ionloss[losses]+addModMass)/charge)

                        valAddKey(all_ions_dict,"z"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+hydrogenAtom)
                        valAddKey(all_ions_dict,"z1"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+(2*hydrogenAtom))
                        valAddKey(all_ions_dict,"z2"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+(3*hydrogenAtom))
                        valAddKey(all_ions_dict,"z3"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+(4*hydrogenAtom))
                    else:
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"x"+losses+"+"+str(charge),"---")

                        valAddKey(all_ions_dict,"z"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"z1"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"z2"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"z3"+losses+"+"+str(charge),"---")
                 
                                

    df = pd.DataFrame(all_ions_dict)
    return df

#This fucntion helps tp accurately idenfy phospho position so that the maximum and minium position can be found to calculate phospho loss
def getPhoPosition(massPosDict):
    phoPos = []
    for key, value in massPosDict.items():
        if "79.96" in str(value):
            phoPos.append(int(key))
    return sorted(phoPos)


def valAddKey(dict1, key, val):
    if key not in dict1.keys():
        dict1[key] = [val]
    else:
        dict1[key].append(val)
    return dict1


def checkAA(pepSeq, check_list):
    update_val = "False"
    for aa in list(pepSeq):
        if aa in check_list:
            update_val = "True"
    return update_val

