|Travis|

Denovonear
----------

This code assesses whether de novo single-nucleotide variants are closer
together within the coding sequence of a gene than expected by chance. We use
local-sequence based mutation rates to account for differential mutability of
regions. The default rates are per-trinucleotide based see `Nature Genetics
46:944–950 <http://www.nature.com/ng/journal/v46/n9/full/ng.3050.html>`_, but
you can use your own rates, or even longer sequence contexts, such as 5-mers or
7-mers.

Install
-------

.. code:: bash
    pip install denovonear


Usage
-----

Analyse *de novo* mutations within python:

.. code:: python

    from denovonear.cluster_test import cluster_de_novos
    
    symbol = 'PPP2R5D'
    de_novos = {'missense': [42975003, 42975003, 42975003, 42975013], 'nonsense': []}
    p_values = cluster_de_novos(symbol, de_novos, iterations=1000000)


Pull out site-specific rates by creating Transcript objects, then get the
rates by consequence at each site

.. code:: python

    from denovonear.ensembl_requester import EnsemblRequest
    from denovonear.load_mutation_rates import load_mutation_rates
    from denovonear.load_gene import construct_gene_object
    from denovonear.site_specific_rates import SiteRates
    
    # convenience object to extract transcript coordinates and sequence from Ensembl
    ensembl = EnsemblRequest(cache_folder='cache', genome_build='grch37')
    transcript = construct_gene_object(ensembl, 'ENST00000346085')
    mut_rates = load_mutation_rates()
    
    rates = SiteRates(transcript, mut_rates)
    
    # rates are stored by consequence, but you can iterate through to find all
    # possible sites in and around the CDS:
    for cq in ['missense', 'nonsense', 'splice_lof', 'synonymous']:
        for site in rates[cq]:
            site['pos'] = transcript.get_position_on_chrom(site['pos'], site['offset'])
    
    # or if you just want the summed rate
    rates['missense'].get_summed_rate()


You can also analyse de novo clustering via the denovonear command:

.. code:: bash

    denovonear cluster \
       --in data/example_de_novos.txt \
       --out output.txt

That command uses a minimal example de novo input file, included in the git
repository. The input is a tab-separated file with a line for each de novo
event. The columns are HGNC symbol, chromosome, position, VEP consequence for
the variant, and whether the de novo is a SNP or indel (the analysis excludes
indels).

Other options are:

* ``--rates PATH_TO_RATES``
* ``--cache-folder PATH_TO_CACHE_DIR``
* ``--genome-build "grch37" or "grch38" (default=grch37)``

The optional rates file is a table separated file with three columns: 'from',
'to', and 'mu_snp'. The 'from' column contains DNA sequence (where the length
is an odd number) with the base to change at the central nucleotide. The 'to'
column contains the sequence with the central base modified. The 'mu_snp' column
contains the probability of the change (as per site per generation).

The cache folder defaults to making a folder named "cache" within the working
directory. The genome build indicates which genome build the coordinates of the
de novo variants are based on, and defaults to GRCh37.

Identify transcripts containing de novo events
----------------------------------------------

You can identify transcripts containing de novos events with the
``identify_transcripts.py`` script. This either identifies all transcripts for a
gene with one or more de novo events, or identifies the minimal set of
transcripts to contain all de novos (where transcripts are prioritised on the
basis of number of de novo events, and length of coding sequence). Transcripts
can be identified with:

.. code:: bash

    denovonear transcripts \
        --de-novos data/example_de_novos.txt \
        --out output.txt \
        --all-transcripts

Other options are:

* ``--minimise-transcripts`` in place of ``--all-transcripts``, to find the minimal
  set of transcripts
* ``--genome-build "grch37" or "grch38" (default=grch37)``

Gene or transcript based mutation rates
---------------------------------------

You can generate mutation rates for either the union of alternative transcripts
for a gene, or for a specific Ensembl transcript ID with the
``construct_mutation_rates.py`` script. Lof and missense mutation rates can be
generated with:

.. code:: bash

    denovonear rates \
        --genes data/example_gene_ids.txt \
        --out output.txt

The tab-separated output file will contain one row per gene/transcript, with
each line containing a transcript ID or gene symbol, a log10 transformed
missense mutation rate, a log10 transformed nonsense mutation rate, and a log10
transformed synonymous mutation rate.

.. |Travis| image:: https://travis-ci.org/jeremymcrae/denovonear.svg?branch=master
    :target: https://travis-ci.org/jeremymcrae/denovonear
