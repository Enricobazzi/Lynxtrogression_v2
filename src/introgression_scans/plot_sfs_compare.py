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
        return 'different'

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

#Data = []
#pop_pairs = ['lpa-wel', 'lpa-sel']
#migs = ['none', 'ab', 'ba']
#
#for pop_pair in pop_pairs:
#    scenarios = get_scenarios(pop_pair)
#    for scenario in scenarios:
#        data_file = f"data/demographic_inference/{pop_pair}_{scenario}.fs"
#        data = dadi.Spectrum.from_file(data_file)
#        for mig in migs:
#            comp_type = get_comp_type(scenario, mig)
#            sim_file = f"data/introgression_scans/simulations_withM/{pop_pair}.{mig}.concatenated.msOut"
#            sim = dadi.Spectrum.from_ms_file(sim_file, average=True, mask_corners=True)
#            sim = sim.fold()
#            # plot_comparison(data, sim, pop_pair, scenario, mig, comp_type)
#            ll = compute_log_likelihood(data, sim)
#            record = create_record(pop_pair, scenario, mig, comp_type, ll)
#            Data.append(record)
#
#df = pd.DataFrame(Data)
#df.to_csv("plots/revision1/compare_sfs/sfs_comparison_log_likelihoods.csv", index=False)
df = pd.read_csv("plots/revision1/compare_sfs/sfs_comparison_log_likelihoods.csv")
# boxplot of log-likelihoods by comparison type
fig, ax = plt.subplots(figsize=(10, 5))
#df.boxplot(column='log_likelihood', by='comparison_type', grid=False, showfliers=True, widths=0.5, boxprops=dict(linewidth=1), medianprops=dict(linewidth=1))
boxprops = dict(linestyle='-', linewidth=1, color='black')
medianprops = dict(linestyle='-', linewidth=1, color='red')
bp = ax.boxplot([df[df['comparison_type'] == 'same direction']['log_likelihood'].values,
                  df[df['comparison_type'] == 'different direction']['log_likelihood'].values],
                 labels=['Same Scenario', 'Different Scenario'],
                 boxprops=boxprops, medianprops=medianprops, patch_artist=False,
                 showfliers=True, widths=0.5)
ax.set_title('Observed Data vs Simulations\n')
ax.set_xlabel('Introgression Scenario')
ax.set_ylabel('Log-Likelihood')
plt.tight_layout()
plt.savefig("plots/revision1/compare_sfs/log_likelihoods_boxplot.pdf")
plt.close()
