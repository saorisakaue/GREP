#!/usr/bin/python
# this will output all drugs and all targets registered in two databases.
import sys
import csv
args = sys.argv #defining args in the command line

# drug bank files

DB = open("./DrugBank/pharmacologically_active.csv", 'rt')
DB_ID = open("./DrugBank/drug_links.csv", 'rt')
UniID = open("HUMAN_9606_idmapping.dat", 'rt')

#DBID_KEGG_Commonname = {}
DBID_Commonname = {}
ID = csv.reader(DB_ID)
header = next(ID)
for row in ID:
	DBID_Commonname[row[0]] = row[1]
DB_ID.close()

Prot_Gene = {}
f = UniID.readlines()
for line in f:
	line = line.rstrip().split("\t")
	if line[1] == "Gene_Name":
		Prot_Gene[line[0]] = line[2] # UnipritID to gene sybbol
UniID.close()

Commonname_target = {}
f = csv.reader(DB)
header = next(f)
for row in f:
	if (row[11] == 'Human') or (row[11] == 'Homo sapiens'):
		DBdrugs =row[12].split("; ")
		for item in DBdrugs:
			if row[2] != "":
				Commonname_target.setdefault(DBID_Commonname[item].lower(),[]).append((row[2],row[0]))
			else:
				if row[0] in Prot_Gene.keys():
					Commonname_target.setdefault(DBID_Commonname[item][0].lower(),[]).append((Prot_Gene[row[0]],row[0]))
DB.close()

# TTD files
DL = open("./TTD/P1-01-TTD_download.txt", 'rt')
Uni = open("./TTD/P2-01-TTD_uniprot_all.txt", 'rt')
Uni_prev = open("./TTD/TTD_uniprot_success.txt", 'rt')

TTD_Uniprot = {}
TName_Uniprot_prev = {}

f = Uni.readlines()
for line in f:
	line = line.rstrip().split("\t")
	uniprot = line[3].split(" ")[0]
	TTD_Uniprot.setdefault(line[0],[]).append(uniprot)
Uni.close()

f = Uni_prev.readlines()
for line in f:
	line = line.rstrip().split("\t")
	targetname = line[2].lower()
	uniprot = line[3]
	TName_Uniprot_prev.setdefault(targetname,[]).append(uniprot)
Uni_prev.close()

TTD_targetname = {}
f = DL.readlines()
for line in f:
	line = line.rstrip().split("\t")
	if line[1] == "Name":
		TTD_targetname[line[0]] = line[2].lower()
	if line[1] == "Drug(s)":
		if line[0] in TTD_Uniprot:
			target_uniprot = TTD_Uniprot[line[0]]
			for prot in target_uniprot:
				if prot in Prot_Gene:
					Commonname_target.setdefault(line[2].lower(), []).append((Prot_Gene[prot],prot))
		else:
			if line[0] in TTD_targetname:
				targetname = TTD_targetname[line[0]]
				if targetname in TName_Uniprot_prev:
					for prot in TName_Uniprot_prev[targetname]:
						if prot in Prot_Gene:
							Commonname_target.setdefault(line[2].lower(), []).append((Prot_Gene[prot],prot))
DL.close()

# remove duplicate
# get all drug-gene relationship
OUT = open("DrugBank_TTD_all_targetgene_information.txt","w")
for drug_name in Commonname_target.keys():
	Commonname_target[drug_name] = list(set(Commonname_target[drug_name]))
	for target in Commonname_target[drug_name]:
		print("{0}\t{1}".format(drug_name,target[0]), file = OUT)
