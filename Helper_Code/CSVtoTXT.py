import csv
from pathlib import Path


# Change these to match your filenames
data_folder = Path('C:\\Users\\brgie\\Documents\\EAE 130')
input_file = data_folder / 'Component_Weights.csv'
output_file = data_folder / 'Ceres100.txt'

with open(input_file, 'r', newline='') as csvfile, open(output_file, 'w') as txtfile:
    reader = csv.reader(csvfile)
    for row in reader:
        txtfile.write(' '.join(row) + '\n')