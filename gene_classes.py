from Bio import pairwise2
from typing import Dict, List


class Gene: 
    def __init__(self, sequence_id:str, sequence_description:str, sequence:str):
        self.__sequence_id = sequence_id
        self.__sequence_description = sequence_description
        self.__sequence = sequence

    @property
    def sequence_id(self):
        return self.__sequence_id
    
    @property
    def sequence_description(self):
        return self.__sequence_description
    
    @property
    def sequence(self):
        return self.__sequence
    
    
    
class MitochondrialDNA(Gene):
     
    def extract_subsequences(self, start=None, end=None):
        """Extract a subsequence from the DNA sequence"""
        sequence = self.sequence
        start = 0 if start is None else start
        end = len(sequence) if end is None else end
        return sequence[start:end]
    
    
    #return the lenght of the frequence
    def sequence_length(self):
        return len(self.sequence)
    
    
    #return all the positions where a certain motif occur
    def find_motif(self, motif: str):

        motif = motif.upper()
        positions = []
        loc = 0
        while True:
            loc = self.sequence.find(motif, loc)        #gives back the location of where it found the motif
            if loc == -1:
                break
            positions.append(loc)
            loc += 1
        return positions
    
    
    #percentage of GC content
    def GC_content_percentage(self):
        sequence=self.sequence
        GC_count=sequence.count('G') + sequence.count('C')
        return f"{(GC_count/len(sequence))*100}%"

    
    #absolute frequency of each nucleotide base
    def base_composition(self):

        return {
            'A': self.sequence.count('A'),
            'C': self.sequence.count('C'),
            'G': self.sequence.count('G'),
            'T': self.sequence.count('T')
        }

    
    #relative frequency of each nucleotide base
    def base_composition_percentage(self):

        return {
            'A': (self.sequence.count('A')/len(self.sequence)*100),
            'C': (self.sequence.count('C')/len(self.sequence)*100),
            'G': (self.sequence.count('G')/len(self.sequence)*100),
            'T': (self.sequence.count('T')/len(self.sequence)*100)
        }
        
        

class BioPythonAligner:
    def align_sequences(self, seq1: str, seq2: str, method: str = 'global') -> Dict:
        """
        Align two sequences using Biopython's pairwise2
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            method: 'global' or 'local'
            
        Returns:
            Dictionary with alignment results
        """
        # Perform alignment
        align_func = pairwise2.align.globalms if method == 'global' else pairwise2.align.localms
        alignments = align_func(seq1, seq2, 2, -1, -0.5, -0.1)  # Default scoring
        
        if not alignments:
            raise ValueError("No alignments found")

        # Get best alignment
        aligned_seq1, aligned_seq2, score, start, end = alignments[0]

        # Calculate matches and gaps
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
        gaps = aligned_seq1.count('-') + aligned_seq2.count('-')
        
        return {
            'aligned_seq1': aligned_seq1,
            'aligned_seq2': aligned_seq2,
            'score': score,
            'matches': matches,
            'gaps': gaps,
            'length': len(aligned_seq1),
            'percent_identity': (matches / len(aligned_seq1)) * 100
        }

    def format_alignment(self, alignment_result: Dict, window_size: int = 60) -> List[Dict]:
        """Format alignment into blocks for display"""
        seq1 = alignment_result['aligned_seq1']
        seq2 = alignment_result['aligned_seq2']
        blocks = []

        for i in range(0, len(seq1), window_size):
            block_seq1 = seq1[i:i+window_size]
            block_seq2 = seq2[i:i+window_size]
            
            # Generate match line for current block
            match_line = ''.join('|' if a == b else ' ' if a == '-' or b == '-' else '.'
                               for a, b in zip(block_seq1, block_seq2))
            
            blocks.append({
                'position': i,
                'seq1': block_seq1,
                'match': match_line,
                'seq2': block_seq2
            })

        return blocks

    
    
