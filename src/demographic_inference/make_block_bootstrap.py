import dadi
import argparse

def parse_args():
    """
    Parse command line arguments: vcf, popmap, out_dir, n_bootstraps
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=str, help="VCF file", required=True)
    parser.add_argument("--popmap", type=str, help="PopMap file", required=True)
    parser.add_argument("--pop_pair", type=str, help="Population Pair: (lpa-wel, lpa-sel, lpa-eel)", required=True)
    parser.add_argument("--out_dir", type=str, help="Output directory", required=True)
    parser.add_argument("--n_bootstraps", type=int, help="Number of bootstraps", default=100)
    parser.add_argument("--chunk_size", type=int, help="Chunk size", default=200_000)
    return parser.parse_args()

def main(vcf, popmap, pop_pair, out_dir, n_bootstraps, chunk_size):
    """
    Create n_bootstraps of the VCF file
    """
    pop_ids = pop_pair.split("-")
    pop_ids.reverse()
    dd = dadi.Misc.make_data_dict_vcf(
        vcf_filename = vcf, popinfo_filename = popmap, filter=True, flanking_info=[None, None]
        )
    chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
    if pop_pair == "lpa-wel":
        ns = [40, 44]
    elif pop_pair == "lpa-eel":
        ns = [38, 44]
    elif pop_pair == "lpa-sel":
        ns = [24, 44]
    boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, n_bootstraps, pop_ids, ns, polarized = False)
    for n in range(len(boots)):
        boots[n].to_file(f'{out_dir}/{pop_pair}_{n}.fs')
    return

if __name__ == "__main__":
    args = parse_args()
    print(f"VCF file: {args.vcf}")
    print(f"PopMap file: {args.popmap}")
    print(f"Populations IDs: {args.pop_pair}")
    print(f"Output directory: {args.out_dir}")
    print(f"Number of bootstraps: {args.n_bootstraps}")
    print(f"Chunk size: {args.chunk_size}")
    main(args.vcf, args.popmap, args.pop_pair, args.out_dir, args.n_bootstraps, args.chunk_size)
