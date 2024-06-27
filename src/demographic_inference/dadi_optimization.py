import gadma2_best_dadi_models as models
import dadi
import argparse
import os

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--boot', type=str, required=False, help='Bootstrapping file from make_block_bootstrap')
    parser.add_argument('--vcf', type=str, required=False, help='VCF file with the data')
    parser.add_argument('--popinfo', type=str, required=False, help='Population info file')
    parser.add_argument('--model', type=str, required=True)
    parser.add_argument('--pop_pair', type=str, required=True)
    parser.add_argument('--csv', type=str, required=True)
    return parser.parse_args()

def optimization_hyperparams():
    # The number of grid points to use in the optimization
    pts = [50, 60, 70]
    # The maximum number of iterations for the optimization
    maxiter = 20
    return pts, maxiter

def get_model(model_name: str):
    try:
        return dadi.Numerics.make_extrap_log_func(getattr(models, f"model_func_{model_name}"))
    except AttributeError:
        raise ValueError(f"Model {model_name} not found in gadma2_best_dadi_models.py")

def get_boot_file(filename: str):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File {filename} not found")
    else:
        return filename
def get_vcf_file(filename: str):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File {filename} not found")
    else:
        return filename

def get_popinfo_file(filename: str):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File {filename} not found")
    else:
        return filename

def get_popids(pop_pair: str):
    popids = pop_pair.split('-')
    popids.reverse()
    return popids

def get_projections(pop_pair: str):
    if pop_pair == 'lpa-wel':
        return [40, 44]
    elif pop_pair == 'lpa-sel':
        return [24, 44]
    elif pop_pair == 'lpa-eel':
        return [38, 44]
    else:
        raise ValueError(f"Projections for population pair {pop_pair} not found")

def get_data_vcf(vcf_file:str, popinfo: str, pop_pair: str):
    vcf_file = get_vcf_file(vcf_file)
    popinfo_file = get_popinfo_file(popinfo)
    pop_ids = get_popids(pop_pair)
    projections = get_projections(pop_pair)
    dd = dadi.Misc.make_data_dict_vcf(vcf_filename=vcf_file, popinfo_filename=popinfo_file, filter=True, flanking_info=[None, None])
    data = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections=projections, polarized=False)
    return data

def get_data_boot(boot:str):
    data = dadi.Spectrum.from_file(boot)
    return data

def get_p0(model_name: str):
    try:
        return getattr(models, f"p0_{model_name}")()
    except AttributeError:
        raise ValueError(f"Model {model_name} not found in gadma2_best_dadi_models.py")

def get_pnames(model_name: str):
    try:
        return getattr(models, f"pnames_{model_name}")()
    except AttributeError:
        raise ValueError(f"Model {model_name} not found in gadma2_best_dadi_models.py")

def main() -> None:
    # read arguments
    args = parse_arguments()
    print(f"Optimizing model {args.model} for population pair {args.pop_pair}\n")
    
    # add log file
    log_file = args.csv
    print(f'Log file: {log_file}')

    # read parameter names
    pnames = get_pnames(args.model)
    pnames_str = ','.join(pnames)

    # write header to log file
    with open(log_file, 'w') as f:
        f.write(f"pop_pair,model,log-likelihood,{pnames_str}\n")
    
    # read data
    print("Reading data...")
    if args.boot:
        data = get_data_boot(args.boot)
    elif args.vcf and args.popinfo:
        data = get_data_vcf(args.vcf, args.popinfo, args.pop_pair)
    else:
        raise ValueError("Either --boot or --vcf and --popinfo must be provided")
    
    # get model function
    model_func = get_model(args.model)
    
    # get optimization hyperparameters
    pts, maxiter = optimization_hyperparams()
    print(f'Optimization hyperparameters: pts={pts}, maxiter={maxiter}')
    
    # get model's initial parameters
    p0 = get_p0(args.model)
    print(f'Initial parameters: {p0}')
    
    # sim model with initial parameters
    model = model_func(p0, data.sample_sizes, pts)
    
    # get initial log-likelihood
    ll_model = dadi.Inference.ll_multinom(model, data)
    print(f'Initial log-likelihood: {ll_model}')
    
    # write initial conditions to log file
    with open(log_file, 'a') as f:
        f.write(f"{args.pop_pair},{args.model},{ll_model},{','.join(str(val) for val in p0)}\n")

    # run optimization
    print(f'Optimization of {args.boot}')
    
    # bound initial parameters
    fixed_p = [p if p == 0.0 else None for p in p0]
    lower_bound = [p * 0.1 for p in p0]
    upper_bound = [p * 10 for p in p0]
    
    # optimize
    opt_params = dadi.Inference.optimize_log(
        p0 = p0, data = data,
        model_func = model_func,
        lower_bound = lower_bound, upper_bound = upper_bound, fixed_params = fixed_p,
        maxiter = maxiter, pts = pts, verbose = 1
        )
    
    # calculate log-likelihood
    model = model_func(opt_params, data.sample_sizes, pts)
    new_ll_model = dadi.Inference.ll_multinom(model, data)
    print(f'optimization {args.boot} log-likelihood: {new_ll_model}')
    print(f'optimization {args.boot} parameters: {opt_params}')
    
    # write results to log file
    with open(log_file, 'a') as f:
        f.write(f"{args.pop_pair},{args.model},{new_ll_model},{','.join(str(val) for val in opt_params.tolist())}\n")
    
    return

## TODO: ADD THETA CALCULATION - NOW IN CONVERSION SCRIPT BUT WOULD BE MORE EFFICIENT HERE

if __name__ == '__main__':
    main()
