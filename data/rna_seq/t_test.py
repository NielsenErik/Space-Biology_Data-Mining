import csv

rows = []
with open("GLDS-48_rna_seq_Normalized_Counts.csv", 'r') as file:
    csvreader = csv.reader(file)
    header = next(csvreader)
    for row in csvreader:
        rows.append(row)
