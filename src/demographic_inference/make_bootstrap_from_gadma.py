"""
Create bootstrap of VCF file using the information contained in the gadma parameter file.

The output is a Frequency Spectrum file in the format of dadi/moments for each bootstrap replicate
saved in the output directory.

Usage:
    make_bootstrap_from_gadma.py <gadma_param_file> <out_dir> <n_bootstraps>

Options:
    <gadma_param_file>  Path to the gadma parameter file.
    <out_dir>           Path to the output directory.
    <n_bootstraps>      Number of bootstraps to create.
"""
import moments
import yaml
import argparse


def parse_gadma_param_file(filename):
    """
    Parse the gadma parameter file to obtain:
    - VCF file name
    - Popinfo file name
    - Populations IDs
    - Projections
    Note: the VCF file and the PopMap are in the same folder as the parameter file.
    """
    with open(filename, "r") as file:
        param = yaml.load(file, Loader=yaml.FullLoader)
    folder = filename.split("/")[:-1]
    folder = "/".join(folder)
    vcf = f'{folder}/{param["Input data"].split(",")[0]}'
    popmap = f'{folder}/{param["Input data"].split(",")[1]}'
    pop_ids = param["Input data"].split(",")[0].split(".")[0].split("-")
    pop_ids.reverse()
    projections = param["Projections"]
    return vcf, popmap, pop_ids, projections

def parse_args():
    """
    Parse command line arguments: gadma_param_file, out_dir, n_bootstraps
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("gadma_param_file", type=str)
    parser.add_argument("out_dir", type=str)
    parser.add_argument("n_bootstraps", type=int)
    return parser.parse_args()

def main(gadma_param_file, out_dir, n_bootstraps):
    """
    Create n_bootstraps of the VCF file using the information contained in the gadma parameter file.
    """
    vcf, popmap, pop_ids, projections = parse_gadma_param_file(gadma_param_file)
    for n in range(0, n_bootstraps):
        data = moments.Misc.make_data_dict_vcf(
            vcf_filename=vcf,
            popinfo_filename=popmap,
            filter=True, flanking_info=[None, None]
        )
        data = moments.Spectrum.from_data_dict(data, pop_ids=pop_ids, projections=projections, polarized=False)
        data.sample().to_file(f"{out_dir}/{'-'.join(pop_ids)}_{n}.fs")
    return

if __name__ == "__main__":
    args = parse_args()
    vcf, popmap, pop_ids, projections = parse_gadma_param_file(args.gadma_param_file)
    print(f"VCF file: {vcf}")
    print(f"PopMap file: {popmap}")
    print(f"Populations IDs: {pop_ids}")
    print(f"Projections: {projections}")
    main(args.gadma_param_file, args.out_dir, args.n_bootstraps)
