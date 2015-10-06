""" functions to load genes, and identify transcripts containing de novos.
"""

from denovonear.transcript import Transcript


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
        dictionary of lengths (in amino acids), indexed by transcript IDs
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
    
    # get the sequence for the identified transcript
    (chrom, start, end, strand, genomic_sequence) = ensembl.get_genomic_seq_for_transcript(transcript_id, expand=10)
    cds_sequence = ensembl.get_cds_seq_for_transcript(transcript_id)
    
    # get the locations of the exons and cds from ensembl
    cds_ranges = ensembl.get_cds_ranges_for_transcript(transcript_id)
    exon_ranges = ensembl.get_exon_ranges_for_transcript(transcript_id)
    
    # start a Transcript object with the locations and sequence
    transcript = Transcript(transcript_id, start, end, strand, chrom, exon_ranges, cds_ranges)
    transcript.add_cds_sequence(cds_sequence)
    transcript.add_genomic_sequence(genomic_sequence, offset=10)
    transcript.fix_coding_sequence_length()
    
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
        try:
            cds_pos = transcript.convert_chr_pos_to_cds_positions(de_novo)
            in_transcript.append(de_novo)
        except ValueError:
            continue
    
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
    transcript_ids = ensembl.get_transcript_ids_for_ensembl_gene_ids(ensembl_genes, [gene_id])
    
    # sometimes we get HGNC symbols that do not match the ensembl rest version
    # that we are currentl using. We can look for earlier HGNC symbols for
    # the gene using the service at rest.genenames.org
    alt_symbols = []
    if len(transcript_ids) == 0:
        alt_symbols = ensembl.get_previous_symbol(gene_id)
        genes = [ensembl.get_genes_for_hgnc_id(symbol) for symbol in alt_symbols]
        genes = [item for sublist in genes for item in sublist]
        ensembl_genes += genes
        symbols = [gene_id] + alt_symbols
        
        transcript_ids = ensembl.get_transcript_ids_for_ensembl_gene_ids(ensembl_genes, symbols)
    
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
        list of Transcript objects for gene, including genomic ranges and sequences
    """
    
    transcripts = minimise_transcripts(ensembl, gene_id, de_novos)
    
    genes = []
    for (transcript_id, total, length) in transcripts:
        gene = construct_gene_object(ensembl, transcript_id)
        genes.append(gene)
    
    if len(genes) == 0:
        raise IndexError("{0}: no suitable transcripts".format(gene_id))
    
    return genes
    
def count_de_novos_per_transcript(ensembl, gene_id, de_novos=[]):
    """ sort out all the necessary sequences and positions for a gene
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        gene_id: HGNC symbol for gene
        de_novos: list of de novo positions, so we can check they all fit in
            the gene transcript
        
    Returns:
        list of (transcript ID, de novo count) tuples, where the de novo count
        shows the number of de novos found in the Ensembl transcript.
    """
    
    transcripts = get_transcript_ids_sorted_by_length(ensembl, gene_id)
    
    # TODO: allow for genes without any coding sequence.
    if len(transcripts) == 0:
        raise IndexError("{0} lacks coding transcripts".format(gene_id))
    
    # count the de novos observed in all transcripts
    counts = []
    for (transcript_id, length) in transcripts:
        try:
            gene = construct_gene_object(ensembl, transcript_id)
            total = len(get_de_novos_in_transcript(gene, de_novos))
            if total > 0:
                counts.append([transcript_id, total, length])
        except ValueError:
            pass
    
    return counts

def minimise_transcripts(ensembl, gene_id, de_novos):
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
        list of [(transcript_id, de novo count, sequence length)] tuples for the
        set of minimal transcripts necessary to contain all de novos.
    """
    
    if len(de_novos) == 0:
        return []
    
    counts = count_de_novos_per_transcript(ensembl, gene_id, de_novos)
    
    if len(counts) == 0:
        return []
    
    # find the transcripts with the most de novos
    max_count = max(item[1] for item in counts)
    transcripts = [item for item in counts if item[1] == max_count]
    
    # find the transcript with the greatest length, should be one transcript
    max_length = max(item[2] for item in transcripts)
    max_transcript = [item for item in transcripts if item[2] == max_length]
    
    # find which de novos occur in the transcript with the most de novos
    gene = construct_gene_object(ensembl, max_transcript[0][0])
    denovos_in_gene = get_de_novos_in_transcript(gene, de_novos)
    
    # trim the de novos to the ones not in the current transcript
    leftovers = [x for x in de_novos if x not in denovos_in_gene]
    
    # and recursively return the transcripts in the current transcript, along
    # with transcripts for the reminaing de novos
    return max_transcript + minimise_transcripts(ensembl, gene_id, leftovers)


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
    
    
