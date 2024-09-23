"""
Create a folded site frequency spectrum file from a VCF file.

Arguments:
--ivcf: Input VCF file
--popinfo: Population information file
--ofs: Output folded site frequency spectrum file
--polarized: True if VCF is polarized, default is False
"""

import dadi
import argparse
from collections import Counter

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ivcf', type=str, required=True,
                        help='Input VCF file')
    parser.add_argument('--popinfo', type=str, required=True,
                        help='Population information file')
    parser.add_argument('--ofs', type=str, required=True,
                        help='Output folded site frequency spectrum file')
    parser.add_argument('--polarized', action='store_true',
                        help='True if VCF is polarized, default is False')
    return parser.parse_args()


def parse_popinfo(popinfo):    
    with open(popinfo, 'r') as f:
        popinfo = f.read().strip().split('\n')
    popinfo = [x.split() for x in popinfo]
    pops = [x[1] for x in popinfo]
    pop_dic = Counter(pops)
    for k in pop_dic.keys():
        pop_dic[k] *= 2
    proj = [pop_dic[x] for x in pop_dic]
    pop_ids = [p for p in pop_dic.keys()]
    return pop_ids, proj

def main():
    args = parse_args()
    pop_ids, proj = parse_popinfo(args.popinfo)
    dd = dadi.Misc.make_data_dict_vcf(vcf_filename = args.ivcf, popinfo_filename = args.popinfo,
                                      filter = True, flanking_info = [None, None])
    data = dadi.Spectrum.from_data_dict(dd, pop_ids = pop_ids, projections = proj, polarized = args.polarized)
    data.to_file(args.ofs)

if __name__ == '__main__':
    main()