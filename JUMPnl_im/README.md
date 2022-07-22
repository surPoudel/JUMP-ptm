# JUMPnl_im #

 * [Introduction](#introduction)
 * [Input Data](#input-data)
 * [Basic Installation](#basic-installation)
 * [JUMPnl_im Commands](#JUMPnl_im-commands)


---

## Introduction ##

JUMPnl_im aims to analyze the precursor neutral loss such as h3po4, hpo3 etc. for phosphorylation or any other neutral loss associated with PTMs. This program also evaluate immonium ions such as acetyllysine for acetylation and other for any PTMs. The main purpose of this program is to validate the JUMPptm results.


[Top of page](#JUMPnl_im)

----


## Input Data ##

JUMPnl_im requires mainly the ms2File and ID.txt file obtained from JUMPptm results and needs to be updated in the parameter file jump_immo_ion_nl_unenrich.params
 - ms2 files - Ideally High Quality spectra; however ms2 file containing all spectra could be used (lower sensitivity) 
 - ID.txt 
 

[Top of page](#JUMPnl_im)

----


## Basic Installation ##

JUMPptm installation should suffice all required packages

----

## JUMPnl_im Commands ##

Once the conda environment (JUMPptm) is activated
1. make a working directory
2. copy the parameter file (jump_immo_ion_nl_unenrich.params) from this directory to the same directory
3. make necessary changes for the parameters (ms2File and idtxt file)
4. Run the command below

```
    python /path of JUMPnl_im/mainProgram.py jump_immo_ion_nl_unenrich.params
```
[Top of page](#JUMPnl_im)


----
Maintainers
----

* To submit bug reports and feature suggestions, please contact

  **Suresh Poudel (suresh.poudel@stjude.org)**

[Top of page](#JUMPptm)

----

