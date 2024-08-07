import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description = "Write fasta files of a particular chromosome from a vcf file")
    parser.add_argument("--ivcf", help = "The input vcf file")
    parser.add_argument("--chr", help = "The chromosome to extract")
    parser.add_argument("--odir", help = "The output directory")
    parser.add_argument("--win_size", help = "The window size")
    parser.add_argument("--step_size", help = "The step size for the sliding window")
    return parser.parse_args()

args = parse_args()

vcf = args.ivcf
chromosome = args.chr
odir = args.odir
window = int(args.win_size)
step = int(args.step_size)

names = []
sequences = {}
positions = []

# parse the vcf file
with open(vcf) as file:
    for line in file:
        if line.startswith("#CHROM"):
            line = line.strip().split()
            for name in line[9:]:
                names.append(f'{name}')
        elif not line.startswith("#"):
            line = line.strip().split()
            if line[0] != chromosome:
                continue
            pos, ref, alt = line[1], line[3], line[4]
            positions.append(pos)
            for i in range(9, len(line)):
                name1 = f'{names[i-9]}.1'
                name2 = f'{names[i-9]}.2'
                if name1 not in sequences:
                    sequences[name1] = []
                if name2 not in sequences:
                    sequences[name2] = []
                if line[i] == "0|0":
                    sequences[name1].append(ref)
                    sequences[name2].append(ref)
                elif line[i] == "0|1":
                    sequences[name1].append(ref)
                    sequences[name2].append(alt)
                elif line[i] == "1|0":
                    sequences[name1].append(alt)
                    sequences[name2].append(ref)
                elif line[i] == "1|1":
                    sequences[name1].append(alt)
                    sequences[name2].append(alt)

# write fasta files
c=0
for i in range(0, len(sequences[name1]), step):
    print(f'writing fasta {c}', end='\r')
    fasta = f'{os.path.abspath(odir)}/{chromosome}_{c}.fa'
    with open(fasta, 'w') as file:
        for name in sequences:
            file.write(f'>{name}:{positions[i:i+window][0]}-{positions[i:i+window][-1]}\n')
            file.write(''.join(sequences[name][i:i+window]) + '\n')
    c += 1

print(f'wrote {c} fasta files of {chromosome} to {odir}')