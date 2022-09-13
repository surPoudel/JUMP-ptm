import collections
import sys
import re
import os
import pandas as pd
import numpy as np
from datetime import datetime
import glob
import argparse
import subprocess
import os.path
import time
import math
from os.path import dirname
import fileinput

def isclose(a, b, tol=1):
  a = float(a)
  b = float(b)
#   massError = abs(a-b)*1e6/a
  massError = (b-a)*1e6/a  #calculates ppm error of the a(theoretical) from b(observed)
  if abs(massError) < tol:
    return True
  else:
    return False



#rushtmt_b22_07
def pepXMLOutFileConversion(pepxmlFile, dtas):
    proton = 1.00727646677
    
    outFileDict, outFileDictPrec = dtaIonDictMap(dtas)
    with fileinput.FileInput(pepxmlFile, inplace=True, backup='.cometOriginal') as f:

        for line in f:
            if "<spectrum_query spectrum=" in line:
                pattern = '<spectrum_query spectrum="(\w+\.\d+\.\d{5}\.\d+)".+precursor_neutral_mass="(\d+\.\d+)"'


                x ="No Pattern"
                allPat=re.match(pattern, line.strip())
                try:
                    x = allPat.group(1)
                    neutralMass = allPat.group(2)
                    MH_mass = float(neutralMass) + proton
                    outfileSplit = x.split(".")

                    for val in range(0,len(outFileDictPrec[x])):
                        if isclose(MH_mass, float(outFileDictPrec[x][val])):
                            break
        #             if x == "HL_human.07462.07462.3":
        #                 print (val,"\t",outFileDict[x][val])
        #                 print (len(outFileDictPrec[x]))
        #                 print (outFileDictPrec[x])
                    #add_value = outFileDict[x][val]
                    add_value = "1"
                    y = outfileSplit[0]+"."+outfileSplit[1]+"."+add_value+"."+outfileSplit[-1]
                    print (line.rstrip().replace(x,y))
                except:
                    print (line.rstrip())
                    
            else:
                print (line.rstrip())

                
def dtaIonDictMap(dtas): #dtas is the dta file
    outFileDict = {} #This is the outfile dictionary that have the ppi information
    outFileDictPrec = {} #This stores the precursor mass to get the accuate ppi information
        
    
    f = open(dtas,"r")
    line=f.readlines()

    for x in range(0, len(line),3):
        dta = line[x]
        mass_ms2 = line[x+1]
        ms2_int = line[x+2]
        dta_info = dta.split(".")
        file = dta_info[0]
        scan = dta_info[1]
        ppi = dta_info[2]
        charge = dta_info[3]
        scan5digit = str('%05d' % int(scan)) #convert to 5 digit for example 1000 scan will be 01000
        scanKey = file+"."+scan5digit+"."+scan5digit+"."+charge #resembles the outfile format of comet pep.xml file

        #if int(charge) != 1: #comet does not search +1 charge for .ms2 so excluding that
        if scanKey not in outFileDict.keys(): #this collects ppi in same order as dtas so these can be used
            outFileDict[scanKey] = [ppi]
        else:
            outFileDict[scanKey].append(ppi)

        neutral_mass = dta.split()[-2] #[M+H]+1

        if scanKey not in outFileDictPrec.keys(): #this collects ppi in same order as dtas so these can be used
            outFileDictPrec[scanKey] = [neutral_mass]
        else:
            outFileDictPrec[scanKey].append(neutral_mass)


    return outFileDict, outFileDictPrec
