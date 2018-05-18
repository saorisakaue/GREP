#!/usr/bin/python
import scipy.stats
import argparse
import string
import sys
# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--genelist','-g',default=None, type=str,
    help='Filename of gene list to analyze.',
	required = True
	)
parser.add_argument('--out','-o',default=None, type=str,
    help='Output filename prefix.',
	required = True
	)
parser.add_argument('--all','-a',default='', type=str,
    help='A list of all genes in the scope of your analysis if available.',
	required = False
	)
parser.add_argument('--test','-t',default=None, type=str,
    help='Test mode : select from ATC or ICD.',
	required = True
	)
parser.add_argument('--emitDrugname','-e',default='False', type=str,
    help='If you want to know drug names target of which overlapped with your genes, set True. [Default = False]',
	required = False
	)
args = parser.parse_args()
IN = open(args.genelist,"r")

target_gene = []
f = IN.readlines()
for line in f:
    line = line.rstrip()
    target_gene.append(line)
target_gene = list(set(target_gene))
IN.close()

all_genes = []
ALL = open("./data/DrugBank_TTD_target_genelist.txt","r")
f = ALL.readlines()
for line in f:
    line = line.rstrip()
    all_genes.append(line)
ALL.close()

if args.all != '':
    all_defined = []
    ALL_File = open(args.all,"r")
    f = ALL_File.readlines()
    for line in f:
        line = line.rstrip()
        all_defined.append(line)
    all_defined = list(set(all_defined))
    all_genes = [gene for gene in all_genes if gene in all_defined]
    ALL_File.close()
all_genes = list(set(all_genes))
target_analysis_genes = [gene for gene in target_gene if gene in all_genes]
target_analysis_genes = list(set(target_analysis_genes))

if args.test == "ATC":
    OUT1 = open(args.out+".ATC.large.txt","a")
    OUT2 = open(args.out+".ATC.detail.txt","a")
    ATC_targetgene = {}
    ATC_drug = {}
    ATCdrug_targetgene = {}
    ATC = open("./data/DrugBank_TTD_targets_by_ATC_v2.txt")
    f = ATC.readlines()
    for line in f:
        line = line.rstrip().split("\t")
        if line[2] != "NA":
            ATC_targetgene.setdefault(line[0], []).append(line[2])
            ATC_drug.setdefault(line[0], []).append(line[1])
            ATCdrug_targetgene.setdefault(line[1], []).append(line[2])
    ATC.close()
    annot_atc = {}
    ANNOT = open("./data/ATC_annotation.txt","r")
    f = ANNOT.readlines()
    for line in f:
        line = line.rstrip().split("\t")
        annot_atc[line[0]] = line[1]
    ANNOT.close()
    largeATC_targetgene = {}
    largeATC_drug ={}
    for atc in ATC_targetgene.keys():
        large_group = atc[0]
        for gene in ATC_targetgene[atc]:
            largeATC_targetgene.setdefault(large_group, []).append(gene)
        for drugname in ATC_drug[atc]:
            largeATC_drug.setdefault(large_group, []).append(drugname)

    # analysis for large group
    if args.emitDrugname == "False":
        print("#Group\tGroupName\tOddsRatio\tFisherExactP",file = OUT1)
    elif args.emitDrugname == "True":
        print("#Group\tGroupName\tOddsRatio\tFisherExactP\tTargetGene:DrugNames",file = OUT1)
    else:
        print("You have to specify True or False[default] for --emitDrugname option. Exit.")
        sys.exit()
    largegroup = list(largeATC_targetgene.keys())
    largegroup.sort()
    for group in largegroup:
        this_group_genes = list(set(largeATC_targetgene[group]))
        joined_genes = [gene for gene in this_group_genes if gene in target_gene]
        joined_genes = list(set(joined_genes))
        classified_genes = [gene for gene in this_group_genes if gene in all_genes]
        classified_genes = list(set(classified_genes))
        a = len(joined_genes)
        b = len(classified_genes) - a
        c = len(target_analysis_genes) - a
        d = len(all_genes) - len(target_analysis_genes) - b
        oddsratio, pvalue = scipy.stats.fisher_exact([[a, b], [c, d]],alternative="greater")
        if args.emitDrugname == "False":
            print("{0}\t{1}\t{2}\t{3}".format(group,annot_atc[group],oddsratio, pvalue),file = OUT1)
        elif args.emitDrugname == "True":
            gene_druglist = []
            for gene in joined_genes:
                joined_drugs = []
                for drug in ATCdrug_targetgene:
                    if gene in ATCdrug_targetgene[drug]:
                        joined_drugs.append(drug)
                joined_drugs = list(set(joined_drugs))
                tmp = ",".join(joined_drugs)
                gene_druglist.append("{0}:{1}".format(gene,tmp))
            drugnames = "\t".join(gene_druglist)
            print("{0}\t{1}\t{2}\t{3}\t{4}".format(group,annot_atc[group],oddsratio, pvalue, drugnames),file = OUT1)
    # analysis for detailed group
    if args.emitDrugname == "False":
        print("#Group\tGroupName\tOddsRatio\tFisherExactP",file = OUT2)
    elif args.emitDrugname == "True":
        print("#Group\tGroupName\tOddsRatio\tFisherExactP\tTargetGene:DrugNames",file = OUT2)
    detailedgroup = list(ATC_targetgene.keys())
    detailedgroup.sort()
    for group in detailedgroup:
        this_group_genes = list(set(ATC_targetgene[group]))
        joined_genes = [gene for gene in this_group_genes if gene in target_gene]
        classified_genes = [gene for gene in this_group_genes if gene in all_genes]
        a = len(joined_genes)
        b = len(classified_genes) - a
        c = len(target_analysis_genes) - a
        d = len(all_genes) - len(target_analysis_genes) - b
        oddsratio, pvalue = scipy.stats.fisher_exact([[a, b], [c, d]],alternative="greater")
        if args.emitDrugname == "False":
            print("{0}\t{1}\t{2}\t{3}".format(group,annot_atc[group],oddsratio, pvalue),file = OUT2)
        elif args.emitDrugname == "True":
            gene_druglist = []
            for gene in joined_genes:
                joined_drugs = []
                for drug in ATCdrug_targetgene:
                    if gene in ATCdrug_targetgene[drug]:
                        joined_drugs.append(drug)
                joined_drugs = list(set(joined_drugs))
                tmp = ",".join(joined_drugs)
                gene_druglist.append("{0}:{1}".format(gene,tmp))
            drugnames = "\t".join(gene_druglist)
            print("{0}\t{1}\t{2}\t{3}\t{4}".format(group,annot_atc[group],oddsratio, pvalue, drugnames),file = OUT2)
    OUT1.close()
    OUT2.close()
elif args.test == "ICD":
    annot_icd = {}
    subcat_icd = {}
    ANNOT = open("./data/icd10_annotation_category.txt")
    f = ANNOT.readlines()
    for line in f:
        line = line.rstrip().split("\t")
        annot_icd[line[0]] = line[1]
        subcat_icd[line[0]] = list(line[2:])
    ANNOT.close()
    annot_ttd = {}
    TTD = open("./data/TTD_crossmatching.txt","r")
    f = TTD.readlines()
    for line in f:
        line = line.rstrip().split("\t")
        if line[2] == "DrugName":
            annot_ttd[line[0]] = line[3].lower()
    TTD.close()
    indicatedICD_drug = {}
    ICD = open("./data/TTD_drug_disease.txt","r")
    f = ICD.readlines()
    for line in f:
        line = line.rstrip().split("\t")
        if line[0] != "TTDDRUGID":
            for drugname in line[0].split("-"):
                drugname = annot_ttd[drugname]
                if len(line) == 5:
                    ICD_10 = line[4].split(", ")
                    for item in ICD_10:
                        if item.find("-") > -1:
                            items = item.split("-")
                            startitem_alp = items[0][0]
                            startitem_num = int(items[0][1:].split(".")[0])
                            if items[1] != "O9A":
                                enditem_alp = items[1][0]
                                enditem_num = int(items[1][1:].split(".")[0])
                                if startitem_alp == enditem_alp:
                                    for i in range(startitem_num,enditem_num+1):
                                        ICDcode = startitem_alp + str('{0:02d}'.format(i))
                                        ICDcode = ICDcode.split(".")[0]
                                        indicatedICD_drug.setdefault(ICDcode, []).append(drugname)
                                else:
                                    start = list(string.ascii_uppercase).index(startitem_alp)
                                    end = list(string.ascii_uppercase).index(enditem_alp)
                                    for i in range(startitem_num,100):
                                        ICDcode = startitem_alp + str('{0:02d}'.format(i))
                                        ICDcode = ICDcode.split(".")[0]
                                        indicatedICD_drug.setdefault(ICDcode, []).append(drugname)
                                    if end - start > 1:
                                        for k in range(start+1,end):
                                            for i in range(0,100):
                                                ICDcode = list(string.ascii_uppercase)[k] + str('{0:02d}'.format(i))
                                                ICDcode = ICDcode.split(".")[0]
                                                indicatedICD_drug.setdefault(ICDcode, []).append(drugname)
                                    for i in range(0,enditem_num+1):
                                        ICDcode = enditem_alp + str('{0:02d}'.format(i))
                                        ICDcode = ICDcode.split(".")[0]
                                        indicatedICD_drug.setdefault(ICDcode, []).append(drugname)
                        else:
                            ICDcode = item.split(".")[0]
                            indicatedICD_drug.setdefault(ICDcode, []).append(drugname)
    ICD.close()
    # get all drugs' tareget information
    alldrug_targetgene = {}
    INFO = open("./data/DrugBank_TTD_all_targetgene_information.txt")
    f = INFO.readlines()
    for line in f:
        line = line.rstrip().split("\t")
        alldrug_targetgene.setdefault(line[0], []).append(line[1])
    INFO.close()
    OUT = open(args.out+".ICD.txt","a")
    if args.emitDrugname == "False":
        print("#Group\tGroupName\tOddsRatio\tFisherExactP",file = OUT)
    elif args.emitDrugname == "True":
        print("#Group\tGroupName\tOddsRatio\tFisherExactP\tTargetGene:DrugNames",file = OUT)
    else:
        print("You have to specify True or False[default] for --emitDrugname option. Exit.")
        sys.exit()
    subcats = list(subcat_icd.keys())
    subcats.sort()
    for cat in subcats:
        cat_name = annot_icd[cat]
        this_group_genes = []
        candidate_drugs = []
        for icd in subcat_icd[cat]:
            if icd in indicatedICD_drug:
                for drug in indicatedICD_drug[icd]:
                    if drug in alldrug_targetgene:
                        candidate_drugs.append(drug)
                        for gene in alldrug_targetgene[drug]:
                            this_group_genes.append(gene)
        this_group_genes = list(set(this_group_genes))
        joined_genes = [gene for gene in this_group_genes if gene in target_gene]
        joined_genes = list(set(joined_genes))
        classified_genes = [gene for gene in this_group_genes if gene in all_genes]
        a = len(joined_genes)
        b = len(classified_genes) - a
        c = len(target_analysis_genes) - a
        d = len(all_genes) - len(target_analysis_genes) - b
        oddsratio, pvalue = scipy.stats.fisher_exact([[a, b], [c, d]],alternative="greater")
        if args.emitDrugname == "False":
            print("{0}\t{1}\t{2}\t{3}".format(cat,cat_name,oddsratio, pvalue),file = OUT)
        elif args.emitDrugname == "True":
            gene_druglist = []
            for gene in joined_genes:
                joined_drugs = []
                for drug in candidate_drugs:
                    if gene in alldrug_targetgene[drug]:
                        joined_drugs.append(drug)
                joined_drugs = list(set(joined_drugs))
                tmp = ",".join(joined_drugs)
                gene_druglist.append("{0}:{1}".format(gene,tmp))
            drugnames = "\t".join(gene_druglist)
            print("{0}\t{1}\t{2}\t{3}\t{4}".format(cat,cat_name,oddsratio, pvalue, drugnames),file = OUT)
    OUT.close()
else:
    print("You have to specify ATC or ICD for --test option. Exit.")
    sys.exit()
