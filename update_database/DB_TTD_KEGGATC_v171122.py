#!/usr/bin/python

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
for drug_name in Commonname_target.keys():
	Commonname_target[drug_name] = list(set(Commonname_target[drug_name]))

# give ATC classification
ATC = open("br08303.keg",'rt')
f = ATC.readlines()

subcat = ["A01","A02","A03","A04","A05","A06","A07","A08","A09","A10","A11","A12","A13","A14","A15","A16","B01","B02","B03","B05","B06","C01","C02","C03","C04","C05","C07","C08","C09","C10","D01","D02","D03","D04","D05","D06","D07","D08","D09","D10","D11","G01","G02","G03","G04","H01","H02","H03","H04","H05","J01","J02","J04","J05","J06","J07","L01","L02","L03","L04","M01","M02","M03","M04","M05","M09","N01","N02","N03","N04","N05","N06","N07","P01","P02","P03","R01","R02","R03","R05","R06","R07","S01","S02","S03","V01","V03","V04","V06","V07","V08","V09","V10"]
# split file into parts
subcat_line_start = []
i = 0
for line in f:
	line = line.rstrip().split()
	if line[0] == "B":
		subcat_line_start.append(i)
	i +=1
ATC.close()

# make subgroup information in KEGG
#KEGGcat_KEGGID = {}
KEGGcat_Commonname = {}

for cat_num in range(len(subcat)-1):
	ATC = open("br08303.keg",'rt')
	f = ATC.readlines()
	cat_name = subcat[cat_num]
	start_num = subcat_line_start[cat_num]
	end_num = subcat_line_start[cat_num+1]-1
	i = 0
	for line in f:
		line = line.rstrip().split()
		if start_num <  i and i < end_num:
			if line[0] == "E":
				druginfo = line[2:]
				if not "and" in line and not "combination" in druginfo and not "Combination" in druginfo and not "combinations" in druginfo and not "Combinations" in druginfo:
					drugname = druginfo[0]
					for k in range(1,len(druginfo)):
						if (druginfo[k].find("[") == -1) and (druginfo[k].find("(") == -1):
							drugname += " " + druginfo[k]
					KEGGcat_Commonname.setdefault(cat_name, []).append(drugname.lower())
			#elif line[0] == "F":
				#KEGGcat_KEGGID.setdefault(cat_name, []).append(line[1])
		i += 1
	ATC.close()

ATC = open("br08303.keg",'rt')
f = ATC.readlines()
cat_num = len(subcat)-1
cat_name = subcat[cat_num]
start_num = subcat_line_start[cat_num]
end_num = 13000
i = 0
for line in f:
	line = line.rstrip().split()
	if start_num <  i and i < end_num:
		if line[0] == "E":
			druginfo = line[2:]
			if not "and" in line and not "combination" in druginfo and not "Combination" in druginfo and not "combinations" in druginfo and not "Combinations" in druginfo:
				drugname = druginfo[0]
				for k in range(1,len(druginfo)):
					if (druginfo[k].find("[") == -1) and (druginfo[k].find("(") == -1):
						drugname += " " + druginfo[k]
				KEGGcat_Commonname.setdefault(cat_name, []).append(drugname.lower())
		#elif line[0] == "F":
		#	KEGGcat_KEGGID.setdefault(cat_name, []).append(line[1])
	i += 1
ATC.close()

OUT = open("DrugBank_TTD_targets_by_ATC_v2.txt","a")
for cat in subcat:
	if cat in KEGGcat_Commonname.keys():
		for drug_name in KEGGcat_Commonname[cat]:
			if drug_name in Commonname_target.keys():
				for info in Commonname_target[drug_name]:
					out = cat + "\t" + drug_name + "\t"
					out += info[0] + "\t" + info[1]
					print(out, file = OUT)
			else:
				out = cat + "\t" + drug_name + "\tNA\tNA"
				print(out, file = OUT)


# output drug information
