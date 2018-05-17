#!/usr/bin/python
import scipy.stats
import argparse
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

if args.test == "ATC":
    OUT1 = open(args.out+".ATC.large.txt","a")
    OUT2 = open(args.out+".ATC.detail.txt","a")
    ATC_targetgene = {}
    ATC = open("./data/DrugBank_TTD_targets_by_ATC_v2.txt")
    f = ATC.readlines()
    for line in f:
        line = line.rstrip().split("\t")
        if line[2] != "NA":
            ATC_targetgene.setdefault(line[0], []).append(line[2])
    ATC.close()
    annot_atc = {}
    ANNOT = open("./data/ATC_annotation.txt","r")
    f = ANNOT.readlines()
    for line in f:
        line = line.rstrip().split("\t")
        annot_atc[line[0]] = line[1]
    ANNOT.close()
    largeATC_targetgene = {}
    for atc in ATC_targetgene.keys():
        large_group = atc[0]
        for gene in ATC_targetgene[atc]:
            largeATC_targetgene.setdefault(large_group, []).append(gene)
    target_analysis_genes = [gene for gene in target_gene if gene in all_genes]
    target_analysis_genes = list(set(target_analysis_genes))
    # analysis for large group
    print("#Group\tGroupName\tOddsRatio\tFisherExactP",file = OUT1)
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
        print("{0}\t{1}\t{2}\t{3}".format(group,annot_atc[group],oddsratio, pvalue),file = OUT1)
    # analysis for detailed group
    print("#Group\tGroupName\tOddsRatio\tFisherExactP",file = OUT2)
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
        print("{0}\t{1}\t{2}\t{3}".format(group,annot_atc[group],oddsratio, pvalue),file = OUT2)
    OUT1.close()
    OUT2.close()
