""" class to obtain and prepare sequence data for analysis
"""

from __future__ import print_function
import os

from get_transcript_sequences import GetTranscriptSequence

class PrepareSequenceData(object):
    """ prepares the data for analysis
    """
    
    def __init__(self, hgnc_mapper):
        
        # self.cds_file = cds_filename
        # self.genomic_file = genomic_sequences_filename
        
        self.hgnc_mapper = hgnc_mapper
        self.ensembl = GetTranscriptSequence()
        
        # self.protein_seqs = self.load_fasta(protein_seqs_filename)
        # self.transcripts = self.identify_transcripts(hgnc_mapper)
        
        # self.prepare_transcript_sequences(self.transcripts)
    
    def get_cds_sequence(self, sequence_id):
        """ returns a fasta dict, indexed by ensembl transcript ID
        """
        
        return self.ensembl.get_cds_seq_for_transcript(sequence_id)
        
        # if not hasattr(self, "cds_seq"):
        #     self.cds_seq = self.load_fasta(self.cds_file)
        
        # return self.cds_seq
    
    def get_genomic_sequence(self, sequence_id, expand):
        """ returns a genome sequence for a transcript
        """
        
        # fasta = self.load_fasta(filename, seq_name=sequence_id)
        
        # return fasta[sequence_id]
        
        return self.ensembl.get_genomic_seq_for_transcript(sequence_id, expand)
    
    def get_gene_start_and_end_for_transcript(self, transcript_id, gene_id):
        """ returns start and end points for a transcript
        """
        
        return self.ensembl.get_gene_start_and_end_for_transcript(transcript_id, gene_id)
    
    def get_exon_ranges_for_transcript(self, transcript_id):
        """ returns a genome sequence for a transcript
        """
        
        return self.ensembl.get_exon_ranges_for_transcript(transcript_id)
    
    def get_cds_ranges_for_transcript(self, transcript_id):
        """ returns a genome sequence for a transcript
        """
        
        return self.ensembl.get_cds_ranges_for_transcript(transcript_id)
    
    def load_fasta(self, filename, trim_unavailable=True, seq_name=None):
        """ creates fasta dict, indexed by ensembl transcript ID
        
        This is so we can find the transcript with longest protein sequence
        for a HGNC ID, which then allows us to obtain and use the CDS of the 
        longest transcript.
        """
        
        fasta_dict = {}
        with open(filename) as f:
            for line in f:
                if line.startswith(">"):
                    transcript_id = line.strip()[1:]
                    if seq_name is None or transcript_id == seq_name:
                        fasta_dict[transcript_id] = ""
                else:
                    if seq_name is None or transcript_id == seq_name:
                        fasta_dict[transcript_id] += line.strip()
                
                if seq_name is not None:
                    if len(fasta_dict) > 0 and transcript_id != seq_name:
                        return fasta_dict
        
        transcript_ids = fasta_dict.keys()
        
        # trim out the unavailable sequences, sometimes we don't want to do
        # this, for example when we want a complete list of all the sequences
        # in a fasta file
        if trim_unavailable:
            for transcript_id in transcript_ids:
                if fasta_dict[transcript_id] == "Sequence unavailable":
                    del fasta_dict[transcript_id]
    
        return fasta_dict
    
    def identify_transcript(self, hgnc_id):
        """ Picks a transcript for a HGNC ID. Currently uses the longest.
        """
        
        transcript_ids = self.hgnc_mapper[hgnc_id]
        transcript_id = self.find_longest_transcript(transcript_ids)
        
        return transcript_id
    
    def find_longest_transcript(self, transcript_ids):
        """ for a given HGNC ID, finds the transcript with the longest CDS
        """
        
        max_length = 0
        max_transcript_id = None
        
        for transcript_id in transcript_ids:
            # get the transcript's protein sequence via the ensembl REST API
            seq = self.ensembl.get_protein_seq_for_transcript(transcript_id)
            
            # ignore transcripts without protein sequence
            if seq == "Sequence unavailable":
                continue
            
            # only swap to using the transcript if it is the longest
            if len(seq) > max_length:
                max_length = len(seq)
                max_transcript_id = transcript_id
        
        return max_transcript_id
    
    def prepare_transcript_sequences(self, hgnc_mapper):
        """ obtains the genomic and cDNA sequences from ensembl for the longest 
        transcript for each of the HGNC IDs.
        """
        
        # load the CDS sequences file, so we can check which sequences have 
        # already been requested,so we don't re-request them (saves time)
        current_seqs = set([])
        if os.path.exists(self.cds_file):
            current_seqs = self.load_fasta(self.cds_file, trim_unavailable=False)
            current_seqs = set(current_seqs.keys())
        
        genomic_seq = open(self.genomic_file, "a")
        cds_seq = open(self.cds_file, "a")
        number_run = 0
        for hgnc_id in hgnc_mapper:
            number_run += 1
            transcript_id = hgnc_mapper[hgnc_id]
            
            # don't re-request the transcript sequence if we have obtained it
            # previously
            if hgnc_id in current_seqs:
                continue
            
            genomic = self.ensembl.get_genomic_seq_for_transcript(transcript_id)
            genomic_seq.write(">{0}\n{1}\n".format(hgnc_id, genomic))
            
            cds = self.ensembl.get_cds_seq_for_transcript(transcript_id)
            cds_seq.write(">{0}\n{1}\n".format(hgnc_id, cds))
            
            # record sequences obtained, and occasionally report progress
            current_seqs.add(hgnc_id)
            if number_run % 100 == 0:
                print(number_run, "files obtained of", len(hgnc_mapper))

        genomic_seq.close()
        cds_seq.close()
    
