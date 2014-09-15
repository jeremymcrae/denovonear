## Denovonear

This code assesses whether recurrent de novo single-nucleotide variants lie 
closer together within the coding sequence of a gene than expected by chance. 
We use mutation rates based on the local sequence context to determine the 
expected likelyhood that specific regions of the gene contain mutations. 
Currently, the local sequence context mutation rates are per-trinucleotide 
mutation ratesb provided by Kaitlin Samocha of the Broad Institute.

Obtain the code from the bitbucket repository with:

`git clone https://github.com/jeremymcrae/denovonear.git`

Then you can shift to the code folder, and set up the code with:
```sh
cd  denovonear/
python setup.py install --user
```

### Analysis
Analyse your de novos with:
```sh
python clustering.py \
   --in data/example_de_novos.txt \
   --rates data/forSanger_1KG_mutation_rate_table.txt \
   --out test.txt
```

That command uses a minimal example de novo input file, included in the git 
repository. The input is a tab-separated file with a line for each de novo 
event. The columns are HGNC symbol, chromosome, position, VEP consequence for 
the variant, and whether the de novo is a SNP or indel (the analysis excludes 
indels). 

Other options are:
 * `--deprecated-genes data/deprecated_ddg2p_hgnc_ids.txt`
 * `--cache-folder PATH_TO_CACHE_DIR`
 * `--genome-build "grch37" or "grch38"`

The deprecated gene ID file is a manually generated file for the genes where 
the code failed to retrieve coordinates from Ensembl, so it's built from the 
gene symbols from analysed files to date, so ongoing analyses might pick up 
additional genes with out of date gene symbols.

The cache folder defaults to making a folder named "cache" within the working 
directory. The genome build indicates which genome build the coordinates of the
de novo variants are based on, and defaults to GRCh37.