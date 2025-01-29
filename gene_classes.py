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

    
    
