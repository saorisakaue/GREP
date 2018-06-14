# GREP: Genome for REPositioning drugs `v1.0.0`
GREP can quantify an enrichment of the user-defined set of genes in the target of clinical indication categories and capture potentially repositionable drugs targeting the gene set.
Both can be run in a few seconds!

## Overview

<div align="center">
<img src="https://raw.githubusercontent.com/saorisakaue/GREP/images/images/forWeb_Overview.png" width=60%>
</div>

The user input is a gene list from any source of genomic studies, and GREP tells you (i) what kind of disease categories are pharmaco-genetically associated with the gene set and (ii) what kind of medications can have a potential for being repositioned to another indication.

## Requirements
`GREP` is a commandline `python` software, and the following modules are required.
- scipy
- argparse
- numpy
- pandas

## Installation
In order to get started with `GREP`, you can clone this repo by the following commands,
```{bash}
git clone https://github.com/saorisakaue/GREP
cd ./GREP
```
Or, you can also use `pip` to install
```{bash}
under dev
```

## Usage
```bash
python grep.py --genelist [list_of_genes] --out [output prefix]  --test [ATC or ICD]
```

`./example/megastroke.genes` can be used as genes identified by megastroke consortium.
`./example/RA_trans.genes` can be used as genes identified in RA meta-analysis by Okada et al.

### Options
`--genelist`, `-g`:  Input your list of genes. Text file with one column. [Required]

`--output`, `-o`:  This will be used as an output prefix. [Required]

`--test`, `-t`:  Choose from 'ATC' or 'ICD' according to the drug categorization. [Required]

`--output-drug-name`, `-d`:  If you want to know drug names target of which overlapped with your genes, set this flag (without any arguments after that).[Optional][Default = False]

`--background`, `-b`: A list of all genes in the scope of your analysis if available. This will be used as background genes.[Optional][Default = None]

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/3.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/3.0/88x31.png" /></a><br />
