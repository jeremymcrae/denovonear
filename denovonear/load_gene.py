""" functions to load genes, and identify transcripts containing de novos.
"""

import asyncio
import logging

from denovonear.transcript import Transcript
from denovonear.gencode import Gencode, Gene

from denovonear.ensembl_requester import (get_protein_seq_for_transcript,
    get_genomic_seq_for_transcript, get_cds_seq_for_transcript,
    get_cds_ranges_for_transcript, get_exon_ranges_for_transcript,
    get_genes_for_hgnc_id, get_transcript_ids_for_ensembl_gene_id,
    get_previous_symbol)

def cds_length(transcript):
    '''get length of CDS
    '''
    return transcript.get_coding_distance(transcript.get_cds_end())['pos'] + 1

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
        de_novos: set of chromosome sequence positions for de novo events
    
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
        within_cds = site['pos'] >= 0 and site['pos'] < cds_length(transcript)
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
       list of transcript IDs for a given HGNC symbol.
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
    
    return transcript_ids

async def load_gene(ensembl, gene_id):
    """ sort out all the necessary sequences and positions for a gene
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        gene_id: HGNC symbol for gene
        
    Returns:
        Gene object, containing info for each transcript
    """
    
    tx_ids = await get_transcript_ids(ensembl, gene_id)
    logging.info(f'found tx IDs for {gene_id}: {tx_ids}')
    
    gene = Gene(gene_id)
    for tx_id in tx_ids:
        logging.info(f'constructing: {tx_id}')
        try:
            tx = await construct_gene_object(ensembl, tx_id)
        except ValueError as e:
            logging.info(f'failed: {tx_id}')
            logging.info(e)
            continue
        gene.add_transcript(tx)
    
    if len(gene.transcripts) == 0:
        logging.info(f"{gene_id} lacks coding transcripts")
    
    return gene

def count_de_novos_per_transcript(transcripts, de_novos=[]):
    """ count de novos in transcripts for a gene.
    
    Args:
        transcripts: list of Transcript objects
        de_novos: list of de novo positions, so we can check they all fit in
            the gene transcript
        
    Returns:
        dictionary of lengths and de novo counts, indexed by transcript IDs.
    """
    
    # count the de novos observed in all transcripts
    counts = {}
    for tx in transcripts:
        tx_id = tx.get_name()
        total = len(get_de_novos_in_transcript(tx, de_novos))
        counts[tx_id] = {'n': total, 'len': cds_length(tx)}
    
    return counts

def minimise_transcripts(transcripts, de_novos):
    """ get a set of minimal transcripts to contain all the de novos.
    
    We identify the minimal number of transcripts to contain all de novos. This
    allows for de novos on mutually exclusive transcripts. The transcripts are
    selected on the basis of containing the most number of de novos, while also
    being the longest possible transcript for the gene.
    
    Args:
        transcripts: list of Transcript objects
        de_novos: list of de novo positions
    
    Returns:
        dictionary of lengths and de novo counts, indexed by transcript ID for
        the set of minimal transcripts necessary to contain all de novos.
    """
    
    if len(de_novos) == 0:
        return {}
    
    counts = count_de_novos_per_transcript(transcripts, de_novos)
    
    # find the transcripts with the most de novos
    max_count = max( val["n"] for key, val in counts.items() )
    max_ids = [ key for key, val in counts.items() if val["n"] == max_count ]
    
    if max_count == 0:
        return {}
    
    # find the transcripts with the greatest length, and the most de novos
    max_length = max( counts[k]["len"] for k in max_ids )
    max_transcripts = {x: counts[x] for x in counts if x in max_ids and counts[x]['len'] == max_length}
    
    # find which de novos occur in the transcript with the most de novos
    best = [x for x in transcripts if x.get_name() == next(iter(max_transcripts))][0]
    denovos_in_gene = get_de_novos_in_transcript(best, de_novos)
    
    # trim the de novos to the ones not in the current transcript
    leftovers = [ x for x in de_novos if x not in denovos_in_gene ]
    
    # recursively return the transcripts in the current transcript, along
    # with transcripts for the reminaing de novos
    update = minimise_transcripts(transcripts, leftovers)
    max_transcripts.update(update)
    
    return max_transcripts
