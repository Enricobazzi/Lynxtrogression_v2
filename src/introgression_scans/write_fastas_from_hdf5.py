import h5py
import argparse
import os
import random

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ifile', required=True)
    parser.add_argument('--odir', required=True)
    parser.add_argument('--nseqs', type=int, default=100)
    return parser.parse_args()

def main():
    args = parse_args()
    ifile = h5py.File(args.ifile, 'r')
    keys = list(ifile.keys())
    random.shuffle(keys)
    c = 0
    for key in keys:
        try:
            int(key)
        except ValueError:
            continue
        for i in range(ifile[key]['x_0'].shape[0]):
            pop1 = ifile[key]['x_0'][i][0]
            pop2 = ifile[key]['x_0'][i][1]
            ofile = f'{os.path.abspath(args.odir)}/fake_{key}_{i}.fa'
            with open(ofile, 'w') as of:
                for j in range(len(pop1)):
                    of.write(f'>pop1_{j}\n{"".join(["A" if str(x) == "0" else "T" for x in pop1[j]])}\n')
                for j in range(len(pop2)):
                    of.write(f'>pop2_{j}\n{"".join(["A" if str(x) == "0" else "T" for x in pop2[j]])}\n')
            c += 1
            print(f'nseqs written: {c}', end='\r')
            if c == args.nseqs:
                break
        if c == args.nseqs:
            break

if __name__ == '__main__':
    main()

            