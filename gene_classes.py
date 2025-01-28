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
    
    
    def sequence_length(self):
        """Calculate the length of the sequence"""
        return len(self.sequence)
    
    
    def GC_content(self):
        """Calculate the GC content percentage"""
        sequence = self.sequence
        GC_count = sequence.count('G') + sequence.count('C')
        percentage = (GC_count / len(sequence)) * 100
        return f"{percentage}%"

    
    
