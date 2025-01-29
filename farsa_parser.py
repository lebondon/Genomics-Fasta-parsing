import pandas as pd
from typing import List, Iterator
from pathlib import Path
from gene_classes import MitochondrialDNA

class FastaParser:
    
    def parse_file(file_path="mtDNA_dataset.txt"):
        records = []
        
        with open(file_path) as file:
            current_id = None
            current_desc = None
            current_seq = []
            
            for line in file:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    if current_id:
                        records.append({
                            'sequence_id': current_id,
                            'sequence_description': current_desc,
                            'sequence': ''.join(current_seq)
                        })
                    
                    # Start new sequence
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

    def create_mitochondrial_objects(df: pd.DataFrame):

        mtdna_list = []
        for sequence_id, row in df.iterrows():
            mtdna = MitochondrialDNA(
                sequence_id=sequence_id,
                sequence_description=row['sequence_description'],
                sequence=row['sequence']
            )
            mtdna_list.append(mtdna)
        return mtdna_list