import argparse
import pandas as pd
"""
Write a table with gene IDs and GO terms from a GFF3 file.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--igff3", type = str, required = True)
    parser.add_argument("--otable", type = str, required = True)
    return parser.parse_args()

args = parse_args()

with open(args.igff3) as file:
    lines = file.readlines()

gene_dict = {}

for line in lines:
    line = line.strip().replace(' ', '_').split()
    if line[2] == 'transcript':
        gene_id = line[8].split(';')[1].split('=')[1]
        gene_dict[gene_id] = []
        go_terms = [l for l in line[8].split(';') if l.startswith('GO:')]
        gene_dict[gene_id] = [go for go in go_terms if go not in gene_dict[gene_id]]

for key in gene_dict.keys():
    gene_dict[key] = ";".join(gene_dict[key])

gene_ids = [k for k in gene_dict.keys()]
go_terms = [v for v in gene_dict.values()]

df = pd.DataFrame({'gene_id': gene_ids, 'go_terms': go_terms})
df.to_csv(args.otable, sep = "\t", index = False)
