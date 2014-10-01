## Denovonear

This code assesses whether recurrent de novo single-nucleotide variants lie 
closer together within the coding sequence of a gene than expected by chance. 
We use mutation rates based on the local sequence context to determine the 
expected likelyhood that specific regions of the gene contain mutations. 
Currently, the local sequence context mutation rates are per-trinucleotide 
mutation rates provided by Kaitlin Samocha of the Broad Institute.

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
   --out output.txt
```

That command uses a minimal example de novo input file, included in the git 
repository. The input is a tab-separated file with a line for each de novo 
event. The columns are HGNC symbol, chromosome, position, VEP consequence for 
the variant, and whether the de novo is a SNP or indel (the analysis excludes 
indels). 

Other options are:
 * `--deprecated-genes data/deprecated_ddg2p_hgnc_ids.txt`
 * `--cache-folder PATH_TO_CACHE_DIR`
 * `--genome-build "grch37" or "grch38" (default=grch37)`

The deprecated gene ID file is a manually generated file for the genes where 
the code failed to retrieve coordinates from Ensembl, so it's built from the 
gene symbols from analysed files to date, so ongoing analyses might pick up 
additional genes with out of date gene symbols.

The cache folder defaults to making a folder named "cache" within the working 
directory. The genome build indicates which genome build the coordinates of the
de novo variants are based on, and defaults to GRCh37.

#### Identify transcripts containing de novo events
You can identify transcripts containing de novos events with the 
`identify_transcripts.py` script. This either identifies all transcripts for a
gene with one or more de novo events, or identifyies the minimal set of 
transcripts to contain all de novos (where transcrips are prioritised on the 
basis of number of de novo events, and length of coding sequence). Transcripts
can be identified with:
```sh
python identify_transcripts.py \
    --de-novos data/example_de_novos.txt \
    --out output.txt \
    --rates data/forSanger_1KG_mutation_rate_table.txt \
    --all-transcripts
```
Other options are:
 * `--minimise-transcripts` in place of `--all-transcripts`, to find the minimal
   set of transcripts
 * `--genome-build "grch37" or "grch38" (default=grch37)`

#### Transcript mutation rates
You can generate mutation rates for an Ensembl transcript ID with the 
`construct_mutation_rates.py` script. Lof and missense mutation rates can be 
generated with:
```sh
python construct_mutation_rates.py \
    --transcripts data/example_transcript_ids.txt \
    --rates data/forSanger_1KG_mutation_rate_table.txt \
    --out output.txt
```

The tab-separated output file will contain one row per transcript, with each
line containing a transcript ID, a log10 transformed missense mutation rate, and
a log10 transformed nonsense mutation rate. 