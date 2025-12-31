import dadi
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import numpy as np

def get_scenarios(pop_pair):
    if pop_pair == 'lpa-wel':
        scenarios = ['wel_to_lpa', 'lpa_to_wel', 'real_data']
    elif pop_pair == 'lpa-sel':
        scenarios = ['sel_to_lpa', 'lpa_to_sel', 'real_data']
    return scenarios

def get_comp_type(scenario, mig):
    if (scenario == 'real_data' and mig == 'none') or \
       (mig == 'ab' and scenario.endswith('_to_lpa')) or \
       (mig == 'ba' and scenario.startswith('lpa_to_')):
        return 'same'
    else:
        return 'opposite'

def plot_comparison(data, sim, pop_pair, scenario, mig, comp_type):
    fig = pylab.figure(1)
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(data=data, model=sim, vmin=10, resid_range=50)
    out_file = f"plots/revision1/compare_sfs/{pop_pair}.{scenario}_data.vs.{mig}_sims.{comp_type}.pdf"
    fig.savefig(out_file)

def compute_log_likelihood(data, sim):
    return dadi.Inference.ll_multinom(sim, data)

def create_record(pop_pair, scenario, mig, comp_type, ll):
    return {
        'pop_pair': pop_pair,
        'scenario': scenario,
        'migration': mig,
        'comparison_type': comp_type,
        'log_likelihood': ll
    }

Data = []
pop_pairs = ['lpa-wel', 'lpa-sel']
migs = ['none', 'ab', 'ba']

for pop_pair in pop_pairs:
    scenarios = get_scenarios(pop_pair)
    for scenario in scenarios:
        data_file = f"data/demographic_inference/{pop_pair}_{scenario}.fs"
        data = dadi.Spectrum.from_file(data_file)
        for mig in migs:
            comp_type = get_comp_type(scenario, mig)
            sim_file = f"data/introgression_scans/simulations_withM/{pop_pair}.{mig}.concatenated.msOut"
            sim = dadi.Spectrum.from_ms_file(sim_file, average=True, mask_corners=True)
            sim = sim.fold()
            plot_comparison(data, sim, pop_pair, scenario, mig, comp_type)
            ll = compute_log_likelihood(data, sim)
            record = create_record(pop_pair, scenario, mig, comp_type, ll)
            Data.append(record)

df = pd.DataFrame(Data)
df.to_csv("plots/revision1/compare_sfs/sfs_comparison_log_likelihoods.csv", index=False)
# boxplot of log-likelihoods by comparison type
plt.figure()
df.boxplot(column='log_likelihood', by='comparison_type')
plt.title('Log-Likelihoods by Comparison Type')
plt.suptitle('')
plt.xlabel('Comparison Type')
plt.ylabel('Log-Likelihood')
plt.savefig("plots/revision1/compare_sfs/log_likelihoods_boxplot.pdf")
plt.close()
