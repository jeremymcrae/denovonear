""" Obtains genomic sequence, transcript IDs, gene, exon and CDS coordinates for
transcripts from Ensembl.
"""

import json
import asyncio

def get_base_url(build):
    assert build in ["grch37", "grch38"], f'unknown build: {build}'
    ver = build + '.' if build == "grch37" else ''
    return f'https://{ver}rest.ensembl.org'

async def get_genes_for_hgnc_id(ensembl, hgnc_symbol, build='grch37'):
    """ obtain the ensembl gene IDs that correspond to a HGNC symbol
    
    e.g. http://rest.ensembl.org/xrefs/symbol/homo_sapiens/KMT2A
    """
    headers = {"content-type": "application/json"}
    ext = f"xrefs/symbol/homo_sapiens/{hgnc_symbol}"
    url = f'{get_base_url(build)}/{ext}'
    resp = await ensembl.get(url, headers)
    return [x['id'] for x in json.loads(resp) if x['type'] == 'gene']

async def get_transcript_ids_for_ensembl_gene_id(ensembl, gene_id, symbols, build='grch37'):
    """ fetch the ensembl transcript IDs for a given ensembl gene ID.
    
    Args:
        gene_id: Ensembl gene ID
        symbols: HGNC symbol(s) for gene. We can have multiple symbols if we
            pull previous symbols out.
    """
    chroms = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", \
         "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", \
          "X", "Y"}
    headers = {"content-type": "application/json"}
    
    transcript_ids = []
    ext = f"overlap/id/{gene_id}?feature=transcript"
    url = f'{get_base_url(build)}/{ext}'
    resp = await ensembl.get(url, headers)
    
    for item in json.loads(resp):
        # ignore non-coding transcripts
        if item["biotype"] not in ["protein_coding", "polymorphic_pseudogene"]:
            continue
        
        # ignore transcripts not on the standard chromosomes
        if item["Parent"] != gene_id or item["seq_region_name"] not in chroms or \
                all([s not in item["external_name"] for s in symbols]):
            continue
        transcript_ids.append(item["id"])
    
    return transcript_ids

async def get_previous_symbol(ensembl, hgnc_symbol):
    """ sometimes we get HGNC symbols that do not match the ensembl rest version
    that we are currently using. We can look for earlier HGNC symbols for
    the gene using the service at rest.genenames.org
    
    Args:
        hgnc_symbol: HGNC symbol for the gene (eg "MLL2")
    
    Returns:
        list of deprecated gene symbols (eg ["KMT2A"])
    """
    headers = {"accept": "application/json"}
    ext = f"fetch/symbol/{hgnc_symbol}"
    url = f'https://rest.genenames.org/{ext}'
    resp = await ensembl.get(url, headers)
    data = json.loads(resp)
    
    docs = data["response"]["docs"]
    # strip out any gene entries that have been invalidated
    docs = [ x for x in docs if x["status"] != "Entry Withdrawn"]
    
    if len(docs) == 0:
        return []
    elif len(docs) > 1:
        raise ValueError(f"{hgnc_symbol} has more than one alternate symbol.")
    elif "prev_symbol" in docs[0]:
        return docs[0]["prev_symbol"]
    return []

async def get_genomic_seq_for_transcript(ensembl, transcript_id, expand, build='grch37'):
    """ obtain the sequence for a transcript from ensembl
    """
    headers = {"content-type": "application/json"}
    ext = f"sequence/id/{transcript_id}?type=genomic;expand_3prime={expand};expand_5prime={expand}"
    url = f'{get_base_url(build)}/{ext}'
    resp = await ensembl.get(url, headers)
    
    gene = json.loads(resp)
    if gene["id"] != transcript_id:
        raise ValueError("ensembl gave the wrong transcript")
    
    _, _, chrom, start, end, strand = gene["desc"].split(":")
    start = int(start) + expand
    end = int(end) - expand
    strand = "+" if strand != '-1' else "-"
    
    return (chrom, start, end, strand, gene["seq"])

async def get_cds_seq_for_transcript(ensembl, transcript_id, build='grch37'):
    """ obtain the sequence for a transcript from ensembl
    """
    headers = {"content-type": "text/plain"}
    ext = f"sequence/id/{transcript_id}?type=cds"
    url = f'{get_base_url(build)}/{ext}'
    resp = await ensembl.get(url, headers)
    return resp.decode('utf8')

async def get_protein_seq_for_transcript(ensembl, transcript_id, build='grch37'):
    """ obtain the sequence for a transcript from ensembl
    """
    headers = {"content-type": "text/plain"}
    ext = f"sequence/id/{transcript_id}?type=protein"
    url = f'{get_base_url(build)}/{ext}'
    resp = await ensembl.get(url, headers)
    return resp.decode('utf8')

async def get_ranges_for_tx(ensembl, transcript_id, feature, build='grch37'):
    """ get coordinates for exons (just cds or full exonic)
    """
    headers = {"content-type": "application/json"}
    ext = f"overlap/id/{transcript_id}?feature={feature}"
    url = f'{get_base_url(build)}/{ext}'
    resp = await ensembl.get(url, headers)
    
    ranges = []
    for exon in json.loads(resp):
        if exon["Parent"] != transcript_id:
            continue
        ranges.append((exon["start"], exon["end"]))
    
    return ranges

async def get_exon_ranges_for_transcript(ensembl, transcript_id, build='grch37'):
    """ obtain the sequence for a transcript from ensembl
    """
    return await get_ranges_for_tx(ensembl, transcript_id, 'exon', build)

async def get_cds_ranges_for_transcript(ensembl, transcript_id, build='grch37'):
    """ obtain the sequence for a transcript from ensembl
    """
    return await get_ranges_for_tx(ensembl, transcript_id, 'cds', build)
