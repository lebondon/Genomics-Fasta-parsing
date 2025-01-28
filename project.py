import pandas as pd 
from gene_classes import Gene

# Function to parse FASTA format file into DataFrame
def build_dataframe_from_file(file_path): 
    data = []  # List to store sequence records

    with open(file_path, "r") as file:
        while True:
            line = file.readline()  # Read line by line
            if not line:  # Exit if end of file
                break 

            first_line = line.strip()  # Remove whitespace
            parts = first_line.split(maxsplit = 1)  # Split into ID and description
            sequence_id = parts[0]
            sequence_description = parts[1] if len(parts) > 1 else ""

            # Read sequence lines until next header or EOF
            sequence = ""
            while True:
                line = file.readline().strip()  # Bug: should be file.readline().strip()
                if not line or line[0] == ">":  # Check for EOF or new sequence header
                    break
                sequence += line

            # Add record to data list
            data.append({"sequence_id": sequence_id,
                         "sequence_description": sequence_description,
                         "sequence": sequence})
            
            # Reposition file pointer if new sequence header found
            if line and line[0] == ">":
                file.seek(file.tell() - len(line) - 1)

    # Create DataFrame and set sequence_id as index
    df = pd.DataFrame(data)
    df.set_index("sequence_id", inplace=True)
    return df

# Function to create Gene objects from DataFrame rows
def create_gene_objects(df):
    gene_objects = []
    for sequence_id, row in df.iterrows():
        gene = Gene(sequence_id, row["sequence_description"], row["sequence"])
        gene_objects.append(gene)
    return gene_objects

# Main execution
file_path = "mtDNA_dataset.txt"
df = build_dataframe_from_file(file_path)
print(df)