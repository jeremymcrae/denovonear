![travis](https://travis-ci.org/jeremymcrae/denovonear.svg?branch=master)

### Denovonear

This code assesses whether de novo single-nucleotide variants are closer
together within the coding sequence of a gene than expected by chance. We use
local-sequence based mutation rates to account for differential mutability of
regions. The default rates are per-trinucleotide based see [Nature Genetics
46:944–950](http://www.nature.com/ng/journal/v46/n9/full/ng.3050.html), but
you can use your own rates, or even longer sequence contexts, such as 5-mers or
7-mers.

### Install
```sh
pip install denovonear
```

### Usage
Analyse *de novo* mutations with the CLI tool:

```sh
denovonear cluster \
   --in data/example.grch38.dnms.txt \
   --gencode data/example.grch38.gtf \
   --fasta data/example.grch38.fa \
   --out output.txt
```

explanation of options:
 - `--in`: path to tab-separated table of de novo mutations. See example table below for columns, or `example.grch38.dnms.txt` in data folder.
 - `--gencode`: path to GENCODE annotations in 
   [GTF format](https://www.ensembl.org/info/website/upload/gff.html) for 
   transcripts and exons e.g. 
   [example release](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz). Can be gzipped, or uncompressed.
 - `--fasta`: path to genome fasta, matching genome build of gencode file

If the --gencode or --fasta options are skipped (e.g. `denovonear cluster --in 
INFILE --out OUTFILE`), gene annotations will be retrieved via an ensembl web 
service. For that, you might need to specify `--genome-build grch38` to ensure
the gene coordinates match your de novo mutation coordinates.

* `--rates PATH_TO_RATES`
* `--cache-folder PATH_TO_CACHE_DIR`
* `--genome-build "grch37" or "grch38" (default=grch37)`

The optional rates file is a table separated file with three columns: 'from',
'to', and 'mu_snp'. The 'from' column contains DNA sequence (where the length
is an odd number) with the base to change at the central nucleotide. The 'to'
column contains the sequence with the central base modified. The 'mu_snp' column
contains the probability of the change (as per site per generation).

The cache folder defaults to making a folder named "cache" within the working
directory. The genome build indicates which genome build the coordinates of the
de novo variants are based on, and defaults to GRCh37.

#### Example de novo table

gene_name | chr | pos | consequence | snp_or_indel
 ---      | --- | --- | ---         |  ---
OR4F5 | chr1 | 69500 | missense_variant | DENOVO-SNP
OR4F5 | chr1 | 69450 | missense_variant | DENOVO-SNP

### Python usage

```py
from denovonear.gencode import Gencode
from denovonear.cluster_test import cluster_de_novos

gencode = Gencode('./data/example.grch38.gtf', './data/example.grch38.fa')
symbol = 'OR4F5'
de_novos = {'missense': [69500, 69450, 69400], 'nonsense': []}
p_values = cluster_de_novos(symbol, de_novos, gencode[symbol], iterations=1000000)
```

Pull out site-specific rates by creating Transcript objects, then get the
rates by consequence at each site

```py
from denovonear.rate_limiter import RateLimiter
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.load_gene import construct_gene_object
from denovonear.site_specific_rates import SiteRates

# extract transcript coordinates and sequence from Ensembl
async with RateLimiter(per_second=15) as ensembl:
    transcript = await construct_gene_object(ensembl, 'ENST00000346085')

mut_rates = load_mutation_rates()
rates = SiteRates(transcript, mut_rates)

# rates are stored by consequence, but you can iterate through to find all
# possible sites in and around the CDS:
for cq in ['missense', 'nonsense', 'splice_lof', 'synonymous']:
    for site in rates[cq]:
        site['pos'] = transcript.get_position_on_chrom(site['pos'], site['offset'])

# or if you just want the summed rate
rates['missense'].get_summed_rate()
```

### Identify transcripts containing de novo events

You can identify transcripts containing de novos events with the
`identify_transcripts.py` script. This either identifies all transcripts for a
gene with one or more de novo events, or identifies the minimal set of
transcripts to contain all de novos (where transcripts are prioritised on the
basis of number of de novo events, and length of coding sequence). Transcripts
can be identified with:

```sh
    denovonear transcripts \
        --de-novos data/example_de_novos.txt \
        --out output.txt \
        --all-transcripts
```

Other options are:

* `--minimise-transcripts` in place of `--all-transcripts`, to find the minimal
  set of transcripts
* `--genome-build "grch37" or "grch38" (default=grch37)`

### Gene or transcript based mutation rates
You can generate mutation rates for either the union of alternative transcripts
for a gene, or for a specific Ensembl transcript ID with the
`construct_mutation_rates.py` script. Lof and missense mutation rates can be
generated with:

```sh
denovonear rates \
    --genes data/example_gene_ids.txt \
    --out output.txt
```

The tab-separated output file will contain one row per gene/transcript, with
each line containing a transcript ID or gene symbol, a log10 transformed
missense mutation rate, a log10 transformed nonsense mutation rate, and a log10
transformed synonymous mutation rate.
