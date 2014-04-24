""" functions to load genes
"""

from src.interval import Interval


def get_deprecated_gene_ids(filename):
    """ gets a dict of the gene IDs used during in DDD datasets that have been 
    deprecated in favour of other gene IDs
    """
    
    deprecated = {}
    with open(filename) as f:
        for line in f:
            line = line.strip().split()
            old = line[0]
            new = line[1]
            deprecated[old] = new
    
    return deprecated

def identify_transcript(ensembl, transcript_ids):
    """ for a given HGNC ID, finds the transcript with the longest CDS
    
    Args:
        ensembl: EnsemblRequest object to request sequences and data 
            from the ensembl REST API
        transcript_ids: list of transcript IDs for a single gene
    
    Returns:
        the transcript ID of the longest protein coding transcript in the list
    """
    
    transcripts = {}
    for transcript_id in transcript_ids:
        # get the transcript's protein sequence via the ensembl REST API
        seq = ensembl.get_protein_seq_for_transcript(transcript_id)
        
        # ignore transcripts without protein sequence
        if seq == "Sequence unavailable":
            continue
        
        transcripts[len(seq)] = transcript_id
    
    return transcripts

def construct_gene_object(ensembl, transcript_id):
    """ creates and Interval object for a gene from ensembl databases
    """
    
    print("loading features and sequence")
    # get the sequence for the identified transcript
    (chrom, start, end, strand, genomic_sequence) = ensembl.get_genomic_seq_for_transcript(transcript_id, expand=10)
    cds_sequence = ensembl.get_cds_seq_for_transcript(transcript_id)
    
    # get the locations of the exons and cds from ensembl
    cds_ranges = ensembl.get_cds_ranges_for_transcript(transcript_id)
    exon_ranges = ensembl.get_exon_ranges_for_transcript(transcript_id)
    
    # start an interval object with the locations and sequence
    transcript = Interval(transcript_id, start, end, strand, chrom, exon_ranges, cds_ranges)
    transcript.add_cds_sequence(cds_sequence)
    transcript.add_genomic_sequence(genomic_sequence, offset=10)
    
    return transcript

def check_denovos_in_gene(transcript, de_novos):
    """ make sure that all the  de novos occur in the loaded gene
    """
    
    for pos in de_novos:
        # convert the de novo positions to cds positions, which raises an error
        # if the position is not in the CDS exons
        try:   
            transcript.convert_chr_pos_to_cds_positions(pos)
        except ValueError:
            return False
    
    return True

def load_gene(ensembl, gene_id, de_novos=[]):
    """ sort out all the necessary sequences and positions for a gene
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        gene_id: HGNC symbol for gene
        de_novos: list of de novo positions, so we can check they all fit in 
            the gene transcript
        
    Returns:
        Interval object for gene, including genomic ranges and sequences
    """
    
    print("loading transcript ID")
    ensembl_genes = ensembl.get_genes_for_hgnc_id(gene_id)
    transcript_ids = ensembl.get_transcript_ids_for_ensembl_gene_ids(ensembl_genes, gene_id)
    transcripts = identify_transcript(ensembl, transcript_ids)
    
    # start with the longest transcript
    lengths = sorted(transcripts)[::-1]
    transcript_id = transcripts[lengths[0]]
    
    # TODO: allow for genes without any coding sequence.
    if transcript_id == {}:
        raise ValueError(gene_id + " lacks coding transcripts")
    
    # create a Interval object using the longest transcript, but if we cannot
    # obtain a vaild sequence or coordinates, or the transcript doesn't contain
    # all the de novo positions, run through alternate transcripts in order of
    # length (allows for CSMD2 variant chr1:34071484 and PHACTR1 chr6:12933929).
    transcript = None
    pos = 0
    while transcript is None or (not check_denovos_in_gene(transcript, de_novos) and pos < (len(transcripts) - 1)):
        try:
            transcript_id = transcripts[lengths[pos]]
            pos += 1
            transcript = construct_gene_object(ensembl, transcript_id)
        except ValueError:
            pass
    
    # raise an IndexError if we can't get a transcript that contains all de 
    # novos. eg ZFN467 with chr7:149462931 and chr7:149461727 which are on
    # mutually exclusive transcripts
    if not check_denovos_in_gene(transcript, de_novos):
        raise IndexError(gene_id + " de novos aren't in CDS sequence")
    
    # make sure we have the gene location for loading the conservation scores
    chrom = transcript.get_chrom()
    start = transcript.get_start()
    end = transcript.get_end()
    
    # print("loading conservation scores")
    # # add in phyloP conservation scores
    # conservation_scores = load_conservation_scores(CONSERVATION_FOLDER, chrom, start, end)
    # try:
    #     transcript.add_conservation_scores(conservation_scores)
    # except ValueError:
    #     pass
    
    return transcript
