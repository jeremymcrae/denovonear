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

def get_transcript_lengths(ensembl, transcript_ids):
    """ finds the protein length for ensembl transcript IDs for a gene
    
    Args:
        ensembl: EnsemblRequest object to request sequences and data 
            from the ensembl REST API
        transcript_ids: list of transcript IDs for a single gene
    
    Returns:
        dictionary of transcript IDs, indexed by their length in amino acids
    """
    
    transcripts = {}
    for transcript_id in transcript_ids:
        # get the transcript's protein sequence via the ensembl REST API
        try:
            seq = ensembl.get_protein_seq_for_transcript(transcript_id)
        except ValueError:
            continue
        
        transcripts[transcript_id] = len(seq)
    
    return transcripts

def construct_gene_object(ensembl, transcript_id):
    """ creates an Interval object for a gene from ensembl databases
    """
    
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

def get_de_novos_in_transcript(transcript, de_novos):
    """ get the de novos within the coding sequence of a transcript
    
    Args:
        transcript: Interval object, which defines the transcript coordinates
        de_novos: list of chromosome sequence positions for de novo events
    
    Returns:
        list of de novo positions found within the transcript
    """
    
    in_transcript = []
    for de_novo in de_novos:
        if transcript.in_coding_region(de_novo):
            in_transcript.append(de_novo)
    
    return in_transcript
    
def get_transcript_ids_sorted_by_length(ensembl, gene_id):
    """ gets transcript IDs for a gene, sorted by coding sequence length
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        gene_id: HGNC symbol for gene
    
    Returns:
        list of (transcript ID, length) tuples for each protein coding 
        transcript for a gene. The transcripts are sorted by transcript length,
        with longest first.
    """
    
    print("loading: {0}".format(gene_id))
    ensembl_genes = ensembl.get_genes_for_hgnc_id(gene_id)
    transcript_ids = ensembl.get_transcript_ids_for_ensembl_gene_ids(ensembl_genes, gene_id)
    transcript_lengths = get_transcript_lengths(ensembl, transcript_ids)
    
    # sort by transcript length
    transcripts = sorted(transcript_lengths.items(), key=lambda x: x[1])
    transcripts = list(reversed(transcripts))
    
    return transcripts

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
    
    transcripts = get_transcript_ids_sorted_by_length(ensembl, gene_id)
    
    # TODO: allow for genes without any coding sequence.
    if len(transcripts) == 0:
        raise ValueError("{0} lacks coding transcripts".format(gene_id))
    
    # create a Interval object using the longest transcript, but if we cannot
    # obtain a valid sequence or coordinates, or the transcript doesn't contain
    # all the de novo positions, run through alternate transcripts in order of
    # length (allows for CSMD2 variant chr1:34071484 and PHACTR1 chr6:12933929).
    for (transcript_id, length) in transcripts:
        try:
            gene = construct_gene_object(ensembl, transcript_id)
            if len(get_de_novos_in_transcript(gene, de_novos)) == len(de_novos):
                # halt the loop, since we've found a transcript with all the de
                # novos
                break  
        except ValueError:
            # this error occurs when the transcript sequence from genomic  
            # sequence according to the gene positions, doesn't match the 
            # transcript sequence obtained from ensembl for the transcript ID.
            pass
    
    # raise an IndexError if we can't get a transcript that contains all de 
    # novos. eg ZFN467 with chr7:149462931 and chr7:149461727 which are on
    # mutually exclusive transcripts
    if len(get_de_novos_in_transcript(gene, de_novos)) == len(de_novos):
        raise IndexError("{0}: de novos aren't in CDS sequence".format(gene_id))
    
    return gene

def load_conservation(transcript, folder):
    """ loads the conservation scores at base level for a gene
    """
    
    # make sure we have the gene location for loading the conservation scores
    chrom = transcript.get_chrom()
    start = transcript.get_start()
    end = transcript.get_end()
    
    print("loading conservation scores")
    # add in phyloP conservation scores
    scores = load_conservation_scores(folder, chrom, start, end)
    try:
        transcript.add_conservation_scores(scores)
    except ValueError:
        pass
    
    return transcript
    
    
