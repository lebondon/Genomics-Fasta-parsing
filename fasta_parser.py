import pandas as pd
from typing import List, Iterator
from pathlib import Path
from gene_classes import *


def detect_sequence_type(df:pd.DataFrame):
    test_sequence=df.iloc[0].sequence
    unique_chars = set(test_sequence.upper())
    if unique_chars.issubset(set("ATCGU")):
        return "DNA"
    else:
        return "protein"

class FastaParser:
    
    def parse_file(file_path):
        records = []
        
        with open(file_path) as file:
            current_id = None
            current_desc = None
            current_seq = []
            
            for line in file:
                line = line.strip()
                if not line:
                    continue
                
                #identifying the header sequence
                if line.startswith('>'):
                    #saving the old data if we had it stored    
                    if current_id:
                        records.append({
                            'sequence_id': current_id,
                            'sequence_description': current_desc,
                            'sequence': ''.join(current_seq)
                        })
                    
                    # Start new sequence
                    # Remove '>' and split into id and description (maximum of 1 split)
                    parts = line[1:].split(maxsplit=1)
                    current_id = parts[0]
                    current_desc = parts[1] if len(parts) > 1 else ''
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_id:
                records.append({
                    'sequence_id': current_id,
                    'sequence_description': current_desc, 
                    'sequence': ''.join(current_seq)
                })
        
        df = pd.DataFrame(records)
        df.set_index('sequence_id', inplace=True)
        return df
    
    def create_sequence_objects(df: pd.DataFrame):

        mtdna_list = []
        if detect_sequence_type(df)=="DNA":
            
            for sequence_id, row in df.iterrows():
                
                mtdna = DNASequence(
                sequence_id=sequence_id,
                sequence_description=row['sequence_description'],
                sequence=row['sequence']
                )
                
                mtdna_list.append(mtdna)
                
        else:
            
            for sequence_id, row in df.iterrows():                
                mtdna = DNASequence(
                sequence_id=sequence_id,
                sequence_description=row['sequence_description'],
                sequence=row['sequence']
                )
                    
                mtdna_list.append(mtdna)
            
                
        return mtdna_list
    
    
    