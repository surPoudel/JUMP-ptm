# JUMPptm #

 * [Introduction](#introduction)
 * [Basic Installation](#basic-installation)
 * [Input Data](#input-data)
 * [Sample Data](#sample-data)
 * [JUMPptm Commands](#jumpptm-commands)
 * [Input and Output Data Organization ](#input-and-output-data-organization)

---

## Introduction ##

JUMPptm aims to identify PTM events from unmatched spectra after conventional peptide analysis of the whole proteome. By default, JUMPptm assumes that the whole proteome data has been analyzed by the JUMP suite, which outputs the identification of peptides/proteins, and the high-quality (HQ) unmatched spectra with de novo tags. Such HQ spectra (input #1) are taken as JUMPptm input, and searched against a list of PTMs specified by the user (input #2) using the multi-stage database search strategy using Comet. The PTM list can be guided by the results from an open search using MSFragger. JUMPptm exports the PTM peptide identification with TMT-based quantification.


NOTE : Please use *ptm_pipeline2.py* script to run the entire pipeline for the standalone version. Similarly, *ptm_pipeline.py* script run entire pipeline for Platform LSF to schedule jobs on the computational cluster. As the list of PTMs gets larger or number of mzXML files increases, is is better to use cluster rather than standalone version. For other platform, users may edit the job submission functions.

Note that JUMPptm only supports 64-bit macOS and Linux.



[Top of page](#JUMPptm)


----
## JUMPptm Publication:
  * The manuscript is submitted and this part will be updated later.
  * If you use JUMPptm as part of a publication, please include this reference.

---

## Input Data ##

JUMPptm requires different input for analyzing the PTMs :
 - ms2 files - Ideally High Quality spectra; however ms2 file containing all spectra could be used (lower sensitivity) 
 - PTMs list (may be derived by using open searches or PTMs of interest) -- this can be updated in the parameter file
 - denovo JUMP derived tags file 
 
 
NOTE: To facilitate subsequent stages of PTM analysis after default peptide analysis, the latest JUMP suite (v1.13.1; https://github.com/JUMPSuite/JUMP) outputs three files after analyzing a whole proteome dataset: i) the accepted unique proteins (.fasta format); ii) HQ unmatched spectra generated by the spectrum QC module (.ms2 format); and iii) de novo amino acid tags for each MS2 spectrum (.tags format). The identified proteins are used to generate a customized database to restrict the search space, and the latter two files are input files for JUMPptm. Note that all the HQ spectra are mass-corrected by the JUMP analysis, allowing narrow mass tolerance (e.g., < 6 ppm) for subsequent PTM identification.

## Sample Data ##
 To evaluate the JUMPptm, we provide a data set composed of 2 fractions [sample_data](./sample_data) directory.

| File name  | Description | Total Demo Spectra and tags |
| ------------- | ------------- | ------------- |
| w001.ms2 | High Quality ms2 file | 792 |
| w010.ms2 | High Quality ms2 file | 654 |
| w001.tags | jump derived denovo tags | 792 |
| w010.tags | jump derived denovo tags | 654 |

----

## Basic Installation ##
A basic install is sufficient for multicore laptops, desktops,
workstations and servers.  For a basic install, use the `bootstrap.sh`
script provided in the repo.  Then, 

1. Place the JUMPptm distribution source in the desired location (call
this `<path to JUMPptm>`)
2. Change your working directory to the top level of the JUMPptm install
and run ` bash bootstrap.sh`

----

  * Obtaining JUMPptm source 
You can obtain the latest version of JUMP from git; simple clone the
git repository:

```
    git clone https://github.com/surPoudel/JUMP-ptm.git
```

in the directory _where you would like JUMPptm to be installed_ (call this directory `<path to JUMPptm>`).  Note
that JUMPptm does not support out-of-place installs; the JUMPptm git
repository _is_ the entire installation.  

----

  * Bootstrapping 
To get dependencies installed, we recommend Conda, and we have
provided a bootstrapping script `bootstrap.sh` that downloads all
dependencies and installs them alongside JUMPptm.  Execute

```
    ./bootstrap.sh
```

and it will create a new directory `JUMPptm` in the current working
directory.  The directory `JUMPptm` will contain the conda
environment to be used by JUMPptm.  JUMPptm will be set up to use the PERL
and python interpreters in that environment.

Once `bootstrap.sh` is finished, activate the conda environment


```
    conda activate $PWD/JUMPptm
```


----

## JUMPptm Commands ##

Once the conda environment (JUMPptm) is activated
1. make a working directory
2. keep all the ms2 files and tags file in the same directory
3. copy the parameter file (ptm_pipeline.params) from [parameterFiles](./parameterFiles) to the same directory
4. make necessary changes for the parameters (including PTM searches stages)
5. Run the command below

```
    python /path of JUMPptm/ptm_pipeline.py ptm_pipeline.params
```

----

## Input and Output Data Organization ##

```bash
.
├── Pipeline_Results_OUTPUT_FOLDER          # Output folder that contains pipeline results (suffixed by Pipeline_Results_)
│   ├── comet.params.new          # comet search parameter file (template)
│   ├── ptm_pipeline.log          # ptm pipeline log file
│   ├── ptm_pipeline.params       # ptm pipeline parameter file is copied inside the results folder for record
│   ├── merge_and_consolidation   # Folder that have results after merging and consolidation of PSMS from each stages 
│   │   ├── ID.txt      # merged IDs from all stages
│   │   ├── jump_fq_merged.params         # automatically generated quantification parameter file
│   │   ├── publications                  # merged filtering results for peptide identification
│   │   │   ├── id_all_pep.txt  # merged peptides identification (all proteins)
│   │   │   └── id_uni_pep.txt  # merged peptides identification (unique proteins)
│   │   ├── quan_HH_tmt10_human_comet     # quantification results folder
│   │   │   └── publications    # folder containing quantification of peptides
│   │   │       ├── id_all_pep_quan.txt    # peptides mapped to all proteins
│   │   │       └── id_uni_pep_quan.txt    # peptides mapped to unique protein
│   │   └── results_table
│   │       └── Pan_PTM_Quan_Table.xlsx  # Pan PTM output excel file
│   ├── Stage_1                                     # Stage_1 search results based on parameter file description
│   │   ├── jump_fc_Stage_1_FDR_1.params  # filtering parameter file for stage 1
│   │   ├── stage_1_comet.params          # search comet parameter file (customized automatically by program based on parameter file)
│   │   ├── sum_Stage_1_FDR_1             # filtering result folder
│   │   │   ├── ID.txt          # identified PSMS (Target only) -- at given FDR 
│   │   │   ├── IDwDecoy.txt    # identified PSMS (Target + Decoy) -- at ven FDR
│   │   │   ├── publications    # folder containing tables file for unique peptide and all peptides
│   │   │   │   ├── id_all_pep_1FDR.txt
│   │   │   │   ├── id_all_pep.txt
│   │   │   │   ├── id_all_prot_1FDR.txt
│   │   │   │   ├── id_all_prot.txt
│   │   │   │   ├── id_uni_pep_1FDR.txt
│   │   │   │   ├── id_uni_pep.txt
│   │   │   │   ├── id_uni_prot_1FDR.txt
│   │   │   │   └── id_uni_prot.txt
│   │   │   └── simplified_report
│   │   │       └── id_uni_prot.txt
│   │   └── w001                         # example fraction name searched by the pipeline --  program makes separate folder for each fraction
│   │       ├── comet.params             # search parameter file is copied
│   │       ├── Results_start_scan_0_end_scan_0_min_tag_len_2        # folder containing tags matched intermediate files
│   │       │   ├── expect_minusLog10.pdf
│   │       │   ├── expect_minusLog10.png
│   │       │   ├── spectrum_tag_count.txt
│   │       │   ├── spectrum_unique_tag_table.txt
│   │       │   ├── tag_qc.params
│   │       │   ├── Total_Tag_matched.pdf
│   │       │   ├── Total_Tag_matched.png
│   │       │   ├── w001_reordered_final.pickle
│   │       │   ├── xcorr.pdf
│   │       │   └── xcorr.png
│   │       ├── search_log.txt          # comet search log
│   │       ├── tag_match.log           # tag match program log file
│   │       ├── tag_qc.params           # tag match parameter file
│   │       ├── w001.1.pep.xml          # pep.xml file output with tag match information
│   │       ├── w001.1.txt              # search file in txt file format
│   │       └── w001.ms2 -> /home/spoudel1/conda_work/test/w001.ms2        # input ms2 softlinked to search fraction folder
│   └── Stage_2
│              .
│              .
│              .
│              .
│
├── ptm_pipeline.params                                     # input parameter file
├── w001.1.tags                                             # input tag file
└── w001.ms2                                                # input ms2 file

```
**Pan_PTM_Quan_Table.xlsx** --- concentanated Pan PTM output file

NOTE: The ID.txt file in merge_and_consolidation folder is modified for the sake of concatenation of different stages. The peptides have Z alphabet appended at the Cterminus that designates the stage. Z = Stage_1; ZZ = Stage_2 etc. The original peptide sequence is also retained. This helps in accurate quantification of peptides using the psms that belongs to specfic stage (so we get stagewise unique psms)

----
