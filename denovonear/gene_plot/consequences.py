"""
Copyright (c) 2015 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

class Consequences(object):
    """ class to identify HGVS-like codes for variants
    
    This class relies on a Transcript object: self.transcript, which defines
    exon positions and genomic sequence to get coding positions, codon sequences
    and translations.
    """
    
    def get_conseqence(self, pos, ref, alt):
        """ get a consequence string for a de novo (e.g. R226Q)
        
        Args:
            pos: nucleotide position of the variant (chromosome coordinates)
            ref: reference allele for the variant
            alt: alternate allele of the variant
        
        Returns:
            string denoting the amino acid change at the given residue (e.g
            R226Q for an arginine subsituted to a glutamine at the 226th
            amino acid).
        """
        
        # make sure the alleles match the strand of the gene
        if self.transcript.strand == "-":
            ref = self.transcript.reverse_complement(ref)
            alt = self.transcript.reverse_complement(alt)
        
        pos = self._adjust_position(pos, ref, alt)
        
        if not self.transcript.in_coding_region(pos):
            return self._get_splice_consequence(pos, ref, alt)
        
        if len(alt) != len(ref):
            return self._get_indel_consequence(pos, ref, alt)
        
        return self._get_missense_consequence(pos, ref, alt)
    
    def _adjust_position(self, pos, ref, alt):
        """ adjust the chromosome position for the matched sequence of alleles.
        
        Sometimes we have a ref allele of TGGGG, and an alt allele of TGGGGGG.
        We adjust the nucleotide position for the section of the alleles that
        match, since otherwise we won't have the position for the amino acid
        that is being altered.
        
        Args:
            pos: nucleotide position of the variant (chromosome coordinates)
            ref: reference allele for the variant
            alt: alternate allele of the variant
        
        Returns:
            position adjusted by the number of bases of the alleles that match.
        """
        
        matched = 0
        try:
            while ref[matched] == alt[matched]:
                matched += 1
        except IndexError:
            pass
        
        if self.transcript.strand == "-":
            pos -= matched
        else:
            pos += matched
        
        return pos
    
    def _get_codon_info(self, pos):
        """ find the codon and supporting information for a transcript position
        
        Args:
            pos: nucleotide position of the variant (chromosome coordinates)
        
        Returns:
            dictionary of codon information, such as the codon sequence, codon
            number, amino acid for the codon, location of the variant within the
            codon.
        """
        
        cds_pos = self.transcript.get_position_in_cds(pos)
        
        if self.transcript.strand == "-":
            cds_pos += 1
        
        codon_number = self.transcript.get_codon_number_for_cds_position(cds_pos)
        intra_codon = self.transcript.get_position_within_codon(cds_pos)
        codon = self.transcript.get_codon_sequence(codon_number)
        initial_aa = self.transcript.translate_codon(codon)
        
        info = {"cds_pos": cds_pos, "codon_number": codon_number + 1,
            "codon": codon, "initial_aa": initial_aa, "intra_codon": intra_codon}
        
        return info
    
    def _get_indel_consequence(self, pos, ref, alt):
        """ figure out the HGVS-like code for indels (frameshifts too) variants.
        
        TODO: This function currently does not provide the full indel
        information. For example, frameshift variants should also have the
        distance in residues to the next stop codon (assuming the frameshift
        variant has introduced a nearby stop codon). Deletions and insertions
        should not just list the initial amino acid that differs, they should
        also include the amino acids that have been removed (or introduced).
        
        Args:
            pos: nucleotide position of the variant (chromosome coordinates)
            ref: reference allele for the variant
            alt: alternate allele of the variant
        
        Returns:
            HGVS-like code for a indel variant e.g. R226fs.
        """
        
        codon = self._get_codon_info(pos)
        
        # the frameshift variants are ones where the difference in length
        # between the ref and alt alleles is not divisible by 3
        if abs(len(ref) - len(alt)) % 3 != 0:
            return "{}{}fs".format(initial_aa, codon["codon_number"])
        elif len(ref) > len(alt):
            return "{}{}del".format(initial_aa, codon["codon_number"])
        elif len(ref) < len(alt):
            return "{}{}ins".format(initial_aa, codon["codon_number"])
    
    def _get_splice_consequence(self, pos, ref, alt):
        """ figure out the HGVS-like code for splice site variants.
        
        First find the closest exon, then the code is based off the distance to
        the nearest end. We also need to take the strand into account. See
        http://www.hgvs.org/mutnomen/recs-DNA.html
        
        Args:
            pos: nucleotide position of the variant (chromosome coordinates)
            ref: reference allele for the variant
            alt: alternate allele of the variant
        
        Returns:
            HGVS-like code for a splice site variant e.g. c.1000-2G>A,
            indicating an exon boundary at position 1000 of the transcript,
            where the variant is 2 bp upstream, and the variant is a G
            substituted with an A.
        """
        
        start, end = self.transcript.find_closest_exon(pos)
        
        distance = min(abs(end - pos), abs(start - pos))
        if pos < start:
            cds_pos = self.transcript.get_position_in_cds(start)
            separator = "-"
        else:
            separator = "+"
            cds_pos = self.transcript.get_position_in_cds(end)
        
        cds_pos += 1
        if self.transcript.strand == "-":
            distance -= 1
            separator = {"-": "+", "+": "-"}[separator]
        
        return "c.{}{}{}{}>{}".format(cds_pos, separator, distance, ref, alt)
    
    def _get_missense_consequence(self, pos, ref, alt):
        """ get HGVS-like code for missense variants
        
        Args:
            pos: nucleotide position of the variant (chromosome coordinates)
            ref: reference allele for the variant
            alt: alternate allele of the variant
        
        Returns:
            HGVS-like code for a missense variant e.g. R226Q for an arginine at
            the 226nd amino acid being substituted with a glutamine.
        """
        
        codon = self._get_codon_info(pos)
        
        mutated_codon = list(codon["codon"])
        mutated_codon[codon["intra_codon"]] = alt
        mutated_codon = "".join(mutated_codon)
        mutated_aa = self.transcript.translate_codon(mutated_codon)
        
        return "{}{}{}".format(codon["initial_aa"], codon["codon_number"], mutated_aa)
