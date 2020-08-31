# -*- coding: utf-8 -*-
#!/usr/bin/python3
# Creation of database IARC with column selected for annovar call
# Date:04/10/19
# Laboratoire d'Hematologie
import sys
import argparse
def create dictionnary_ADN():
	ADN = {}
	ADN['A'] = 'T'
	ADN['T'] = 'A'
	ADN['C'] = 'G'
	ADN['G'] = 'C'
	ADN['N'] = "N"

	return ADN

# Function
def create_database(file_in,dico,out):
	out_file = open(out,"w")
	
	colonne_chr = "Chr"
	
	colonne_end = "End"
	sep = "\t"
	with open(file_in,"r") as file_database:
		for ligne in file_database:
			info = ligne.split("\t")
			# ligne d'annotation
			if info[0].startswith("Mut") == False:
				# Recuperation des infos
				start = info[3]
				r = str(info[11])
				a = str(info[12])
				ref = r.replace(" ","",)
				alt = a.replace(" ","",)
				# Transactivation class
				annotation_trans = info[33]
				# Structural class
				annotation_func = info[36]
				# 
				#MUT_ID
				#MUT_ID 	Unique identifier of each gene variation reported in the database. This identifier is used in all datasets (somatic, polymorhisms, germline). 
				ID = info[1]
				Coding_Genomic = info[9]
				if len(ref) == 1 and len(alt) == 1:
					end = int(start)
				# Deletion
				if alt == "NA" :
					alt ="."
				
				# Translation try

				if ref in dico:
					ref = dico[ref]
				if alt in ADN:
					alt = dico[alt]


				ligne_info = "17" + sep + start + sep + str(end) + sep + ref + sep + alt + sep + Coding_Genomic + ";" + annotation_trans + "," + annotation_func
				out_file.write(ligne_info + "\n")
			# Info header
			else:
				colonne_ID = info[1]
				c_genomic = "IARC"
				colonne_annotation = "IARC" + info[33]
				colonne_start = info[3]
				c_ref = info[11]
				c_alt = info[12]
				header_out = "#" + colonne_chr + sep + colonne_start + sep + colonne_end + sep + c_ref + sep + c_alt + sep + "IARC"
				print(header_out)
				# En tete verification
	out_file.close()
			

# *****************************************************
# ********************* Main **************************
# *****************************************************
if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog="Treatment_of_Annotation_v2")
	# Create sub_parsers (one for Treatment_Annotation, other for Combine)
	# File out
	parser.add_argument("-i", "--input",
		dest="primary_database",
		required=True,
		type=str,
		help="directory of result annotation for annovar"
	)
	# File in
	parser.add_argument("-o", "--output",
		dest="transform_database",
		required=True,
		type=str,
		help="File of result annotation for annovar"
	)
	args = parser.parse_args()
	file_in = "/media/t-chu-027/DATAPART1/Database/IARC-TP53/IARC_TP53.csv"
	out_file = "/media/t-chu-027/DATAPART1/Database/humandb_annovar/hg19_IARC.txt"
	dico = create_dictionnary()
	create_database(args.primary_database,dico,args.transform_database)
