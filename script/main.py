#!/usr/bin/env python
import argparse
import os
import scipy.stats
import numpy as np
import pandas as pd

# consts
BASEDIR = os.path.dirname(__file__)
DATADIR = os.path.normpath(os.path.join(BASEDIR, '..', 'data'))
ALL_GENES_FNAME = os.path.join(DATADIR, 'DrugBank_TTD_target_genelist.txt')
ATC_TARGETS_FNAME = os.path.join(DATADIR, 'DrugBank_TTD_targets_by_ATC_v2.txt')
ATC_ANNOT_FNAME = os.path.join(DATADIR, 'ATC_annotation.txt')
ICD_TARGETS_FNAME = os.path.join(DATADIR, 'DrugBank_TTD_targets_by_ICD10.txt')
ICD_ANNOT_FNAME = os.path.join(DATADIR, 'icd10_annotation_category.txt')

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--genelist', '-g', default=None, type=str,
    help='Filename of gene list to analyze.',
    required=True)
parser.add_argument('--out', '-o', default=None, type=str,
    help='Output filename prefix.',
    required=True)
parser.add_argument('--all', '-a', default=None, type=str,
    help='A list of all genes in the scope of your analysis if available.',
    required=False)
parser.add_argument('--test', '-t', default=None, type=str, choices=['ATC', 'ICD'],
    help='Test mode : select from ATC or ICD.',
    required=True)
parser.add_argument('--emitDrugname', '-e', default=False, action='store_true',
    help='If you want to know drug names target of which overlapped with your genes, set this flag.',
    required=False)


# fisher exact test function
def grep(target, key, annot, emitDrugname=False):
    group = target[key].iloc[0]
    this_group_genes = target.TargetGene.unique()
    joined_genes = np.intersect1d(target_gene, this_group_genes)
    classified_genes = np.intersect1d(all_genes, this_group_genes)
    a = len(joined_genes)
    b = len(classified_genes) - a
    c = len(target_analysis_genes) - a
    d = len(all_genes) - len(target_analysis_genes) - b
    oddsratio, pvalue = scipy.stats.fisher_exact([[a, b], [c, d]], alternative="greater")

    out_cols = ['#Group', 'GroupName', 'OddsRatio', 'FisherExactP']
    if emitDrugname:
        out_cols.append('TargetGene:DrugNames')
        gene_druglist = target[target.TargetGene.isin(joined_genes)].groupby('TargetGene').apply(
            lambda x: x.TargetGene.iloc[0] + ':' + ','.join(x.Drug.unique()))
        drugnames = ";".join(gene_druglist)
        ret = [group, annot.Annot[group], oddsratio, pvalue, drugnames]
    else:
        ret = [group, annot.Annot[group], oddsratio, pvalue]
    return pd.Series(ret, out_cols)


if __name__ == '__main__':
    args = parser.parse_args()

    target_gene = pd.read_csv(args.genelist, header=None)[0].unique()
    all_genes = pd.read_table(ALL_GENES_FNAME, header=None)[0]

    if args.all is not None:
        all_defined = pd.read_csv(args.all, header=None)[0].unique()
        all_genes = np.intersect1d(all_genes, all_defined)
    target_analysis_genes = np.intersect1d(all_genes, target_gene)

    if args.test == "ATC":
        ATC = pd.read_csv(ATC_TARGETS_FNAME, header=None, sep='\t',
                names=['Code', 'Drug', 'TargetGene', 'XXX']).dropna(subset=['TargetGene'])
        ATC['large'] = ATC.Code.str.slice(0, 1)
        ATC_ANNOT = pd.read_csv(ATC_ANNOT_FNAME, header=None, sep='\t', index_col=0, names=['Code', 'Annot'])

        # analysis for large group
        ATC.groupby('large').apply(
            grep, annot=ATC_ANNOT, key='large',
            emitDrugname=args.emitDrugname).to_csv(
                args.out + ".ATC.large.txt", index=False, sep='\t')

        # analysis for detailed group
        ATC.groupby('Code').apply(
            grep, annot=ATC_ANNOT, key='Code',
            emitDrugname=args.emitDrugname).to_csv(
                args.out + ".ATC.detail.txt", index=False, sep='\t')

    elif args.test == "ICD":
        annot_icd = {}
        subcat_icd = {}
        with open(ICD_ANNOT_FNAME, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                annot_icd[line[0]] = line[1]
                [subcat_icd.update({v: line[0]}) for v in line[2:]]
        ICD_ANNOT = pd.DataFrame.from_dict(
            annot_icd, orient='index').rename(columns={0: 'Annot'})

        ICD = pd.read_csv(ICD_TARGETS_FNAME, header=None, sep='\t', names=['Code', 'Drug', 'TargetGene'])
        ICD['large'] = ICD.Code.apply(
            lambda x: subcat_icd.setdefault(x, np.nan))
        ICD = ICD.dropna()

        ICD.groupby('large').apply(
            grep, annot=ICD_ANNOT, key='large',
            emitDrugname=args.emitDrugname).to_csv(
                args.out + ".ICD.txt", index=False, sep='\t')
