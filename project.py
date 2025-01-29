import pandas as pd 
from gene_classes import Gene


# Main execution
file_path = "mtDNA_dataset.txt"
df = build_dataframe_from_file(file_path)
print(df)