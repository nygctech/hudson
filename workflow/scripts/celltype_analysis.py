import csv
#import os
import pandas as pd
from pathlib import Path


# Find cell type csv by section
table_path = Path(snakemake.input[0])

# Create dictionary where keys are sections and values are a dictionary of cell types & their fractional counts
# sections = [f.name[0:-4] for f in os.scandir(table_path) if f.name[-4:] == '.csv']

cell_type_counts = {}
cell_type_counts_frac = {}


# Create a dictionary to count category appearances
category_count_dict = {}

# Open the CSV file and read its content
with open(table_path, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)

    # Skip the header row if it exists
    next(csv_reader, None)

    # Iterate through the rows in the CSV file
    for row in csv_reader:
        category = row[1]  # Assuming the category (cell type) is in the second column (index 1)

        if category in category_count_dict:
            category_count_dict[category] += 1
        else:
            category_count_dict[category] = 1                

cell_type_counts = category_count_dict
    
df_celltypecount = pd.DataFrame.from_dict(cell_type_counts, orient='index')
df_celltypecount.reset_index(inplace=True)
df_celltypecount.columns = ['Cell Type', 'Count']
total_count = df['Count'].sum()

# Add column with fractional count
df_celltypecount['Fractional Count'] = df_celltypecount['Count'] / total_count


df_celltypecount.to_csv(snakemake.output[0], index=False)