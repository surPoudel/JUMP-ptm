from Bio import SeqIO
import sys
import glob
import pandas as pd
import os, sys


def make_targetIDs_2_decoy(accepted_protein_list):
	decoy_accepted_protein_list = []
	for prot in accepted_protein_list:
	    decoy_accepted_protein_list.append(prot)
	    decoy_accepted_protein_list.append("##Decoy__"+prot)

	return decoy_accepted_protein_list


def filtered_IDs_to_fasta(database_name, db_folder, accepted_protein_list):

	decoy_accepted_protein_list = make_targetIDs_2_decoy(accepted_protein_list)

	result_file = db_folder+'/custom_filtered_Proteins.fasta'
	fasta_sequences = SeqIO.parse(open(database_name),'fasta')
	with open(result_file, "w") as f:
	  for seq in fasta_sequences:
	    if seq.id in decoy_accepted_protein_list:
	      SeqIO.write([seq], f, "fasta")

	return result_file
