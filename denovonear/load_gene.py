""" functions to load genes, and identify transcripts containing de novos.
"""

import asyncio

from denovonear.transcript import Transcript

from denovonear.ensembl_requester import (get_protein_seq_for_transcript,
    get_genomic_seq_for_transcript, get_cds_seq_for_transcript,
    get_cds_ranges_for_transcript, get_exon_ranges_for_transcript,
    get_genes_for_hgnc_id, get_transcript_ids_for_ensembl_gene_id,
    get_previous_symbol)

async def get_transcript_lengths(ensembl, transcript_ids):
    """ finds the protein length for ensembl transcript IDs for a gene
    
    Args:
        ensembl: EnsemblRequest object to request sequences and data
            from the ensembl REST API
        transcript_ids: list of transcript IDs for a single gene
    
    Returns:
        dictionary of lengths (in amino acids), indexed by transcript IDs
    """
    tasks = [get_protein_seq_for_transcript(ensembl, x) for x in transcript_ids]
    seqs = await asyncio.gather(*tasks, return_exceptions=True)
    
    # associate tx_ids to the protein sequences, but remove failures
    transcripts = dict(zip(transcript_ids, seqs))
    for key in list(transcripts.keys()):
        if not isinstance(transcripts[key], str):
            del transcripts[key]
        else:
            transcripts[key] = len(transcripts[key])
    
    return transcripts

async def construct_gene_object(ensembl, transcript_id):
    """ creates an Transcript object for a gene from ensembl databases
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        transcript_id: string for an Ensembl transcript ID
    
    Returns:
        a Transcript object, containing transcript coordinates and gene and
        transcript sequence.
    
    Raises:
        ValueError if CDS from genomic sequence given gene coordinates and CDS
        retrieved from Ensembl do not match.
    """
    tasks = [get_genomic_seq_for_transcript(ensembl, transcript_id, expand=10),
        get_cds_seq_for_transcript(ensembl, transcript_id),
        get_cds_ranges_for_transcript(ensembl, transcript_id),
        get_exon_ranges_for_transcript(ensembl, transcript_id)]
    (chrom, start, end, strand, genomic), cds_seq, cds, exons = await asyncio.gather(*tasks)
    
    # start a Transcript object with the locations and sequence
    transcript = Transcript(transcript_id, chrom, start, end, strand)
    transcript.set_exons(exons, cds)
    transcript.set_cds(cds)
    transcript.add_cds_sequence(cds_seq)
    transcript.add_genomic_sequence(genomic, offset=10)
    
    return transcript

def get_de_novos_in_transcript(transcript, de_novos):
    """ get the de novos within the coding sequence of a transcript
    
    Args:
        transcript: Transcript object, which defines the transcript coordinates
        de_novos: list of chromosome sequence positions for de novo events
    
    Returns:
        list of de novo positions found within the transcript
    """
    
    in_transcript = []
    for de_novo in de_novos:
        # we check if the de novo is within the transcript by converting the
        # chromosomal position to a CDS-based position. Variants outside the CDS
        # will raise an error, which we catch and pass on. It's better to do
        # this, rather than use the function in_coding_region(), since that
        # function does not allow for splice site variants.
        site = transcript.get_coding_distance(de_novo)
        cds_length = transcript.get_coding_distance(transcript.get_cds_end())
        within_cds = site['pos'] >= 0 and site['pos'] < cds_length['pos']
        if within_cds and (transcript.in_coding_region(de_novo) or abs(site['offset']) < 9):
            in_transcript.append(de_novo)
    
    return in_transcript

def flatten(values):
    return [x for sub in values for x in sub]

async def get_transcript_ids(ensembl, symbol):
    """ gets transcript IDs for a gene.
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        symbol: HGNC symbol for gene
    
    Returns:
        dictionary of transcript ID: transcript lengths for all transcripts
        for a given HGNC symbol.
    """
    
    genes = await get_genes_for_hgnc_id(ensembl, symbol)
    tasks = [get_transcript_ids_for_ensembl_gene_id(ensembl, x, symbol) for x in genes]
    transcript_ids = flatten(await asyncio.gather(*tasks))
    
    # sometimes we get HGNC symbols that do not match the ensembl rest version
    # that we are currentl using. We can look for earlier HGNC symbols for
    # the gene using the service at rest.genenames.org
    symbols = []
    if len(transcript_ids) == 0:
        symbols = await get_previous_symbol(ensembl, symbol)
        tasks = [get_genes_for_hgnc_id(ensembl, s) for s in symbols]
        genes += flatten(await asyncio.gather(*tasks))
        symbols = [symbol] + symbols
        
        tasks = [get_transcript_ids_for_ensembl_gene_id(ensembl, x, symbols) for x in genes]
        transcript_ids = flatten(await asyncio.gather(*tasks))
    
    return await get_transcript_lengths(ensembl, transcript_ids)

async def load_gene(ensembl, gene_id, de_novos=[]):
    """ sort out all the necessary sequences and positions for a gene
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        gene_id: HGNC symbol for gene
        de_novos: list of de novo positions, so we can check they all fit in
            the gene transcript
        
    Returns:
        list of Transcript objects for gene, including genomic ranges and sequences
    """
    
    transcripts = await minimise_transcripts(ensembl, gene_id, de_novos)
    
    genes = []
    for transcript_id in transcripts:
        gene = await construct_gene_object(ensembl, transcript_id)
        genes.append(gene)
    
    if len(genes) == 0:
        raise IndexError("{0}: no suitable transcripts".format(gene_id))
    
    return genes
    
async def count_de_novos_per_transcript(ensembl, gene_id, de_novos=[]):
    """ count de novos in transcripts for a gene.
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        gene_id: HGNC symbol for gene
        de_novos: list of de novo positions, so we can check they all fit in
            the gene transcript
        
    Returns:
        dictionary of lengths and de novo counts, indexed by transcript IDs.
    """
    
    transcripts = await get_transcript_ids(ensembl, gene_id)
    
    # TODO: allow for genes without any coding sequence.
    if len(transcripts) == 0:
        raise IndexError("{0} lacks coding transcripts".format(gene_id))
    
    # count the de novos observed in all transcripts
    counts = {}
    for key in transcripts:
        try:
            gene = await construct_gene_object(ensembl, key)
            total = len(get_de_novos_in_transcript(gene, de_novos))
            if total > 0:
                counts[key] = {}
                counts[key]["n"] = total
                counts[key]["len"] = transcripts[key]
        except ValueError:
            pass
    
    return counts

async def minimise_transcripts(ensembl, gene_id, de_novos):
    """ get a set of minimal transcripts to contain all the de novos.
    
    We identify the minimal number of transcripts to contain all de novos. This
    allows for de novos on mutually exclusive transcripts. The transcripts are
    selected on the basis of containing the most number of de novos, while also
    being the longest possible transcript for the gene.
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        gene_id: HGNC symbol for gene
        de_novos: set of de novo positions
    
    Returns:
        dictionary of lengths and de novo counts, indexed by transcript ID for
        the set of minimal transcripts necessary to contain all de novos.
    """
    
    if len(de_novos) == 0:
        return {}
    
    counts = await count_de_novos_per_transcript(ensembl, gene_id, de_novos)
    
    if len(counts) == 0:
        return {}
    
    # find the transcripts with the most de novos
    max_count = max( val["n"] for key, val in counts.items() )
    transcripts = [ key for key, val in counts.items() if val["n"] == max_count ]
    
    # find the transcripts with the greatest length, and the most de novos
    max_length = max( counts[key]["len"] for key in transcripts )
    tx_ids = [ key for key in transcripts if counts[key]["len"] == max_length ]
    max_transcripts = {x: counts[x] for x in counts if x in tx_ids}
    
    # find which de novos occur in the transcript with the most de novos
    gene = await construct_gene_object(ensembl, next(iter(max_transcripts)))
    denovos_in_gene = get_de_novos_in_transcript(gene, de_novos)
    
    # trim the de novos to the ones not in the current transcript
    leftovers = [ x for x in de_novos if x not in denovos_in_gene ]
    
    # and recursively return the transcripts in the current transcript, along
    # with transcripts for the reminaing de novos
    update = await minimise_transcripts(ensembl, gene_id, leftovers)
    max_transcripts.update(update)
    
    return max_transcripts
