*** data directory is under develoment ***

# GREP: Genome for REPositioning drugs
You can use any gene list to test whether those genes are enriched in certain drug categories.
You are able to know potential drugs that target genes you want to examine.
Both can be run in a few seconds!

## Requirements
- python
- scipy
- argparse
- numpy
- pandas

## Installation
```{bash}
git clone https://github.com/saorisakaue/GREP
cd ./GREP
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

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/3.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/3.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/3.0/">Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Unported License</a>.
