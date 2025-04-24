import demes
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import LogLocator

def get_size_functions(graph):
    """
    Get the size functions from the graph
    """
    pop1_sfs = []
    pop2_sfs = []
    for n in list(reversed(range(3))):
        pop1_sfs.append(graph.demes[2].epochs[n].size_function)
        pop2_sfs.append(graph.demes[3].epochs[n].size_function)
    return pop1_sfs, pop2_sfs

def get_100_points(t0, t1, n0, n1, sf):
    """
    Get 100 points between t0 and t1 for the size function
    """
    t = np.linspace(t0, t1, 100)
    if sf == 'exponential':
        # growth rate
        k = np.log(n1 / n0) / (t1 - t0)
        s = np.array(n0 * np.exp(k * t))
        s = np.array(n0 * np.exp(k * (t - t0)))
    elif sf == 'constant':
        s = np.array([n0 for _ in range(100)])
    return t, s

def get_size_df_row(row, boot_num, graph):
    """
    Get the size DataFrame for a given row
    """
    pop1_sfs, pop2_sfs = get_size_functions(graph)
    e1_t0 = 0
    e1_t1 = row['t4']
    e1_n1_0 = row['nu41']
    e1_n1_1 = row['nu31']
    e1_sf1 = pop1_sfs[0]
    time, size1 = get_100_points(e1_t0, e1_t1, e1_n1_0, e1_n1_1, e1_sf1)
    e1_n2_0 = row['nu42']
    e1_n2_1 = row['nu32']
    e1_sf2 = pop2_sfs[0]
    size2 = get_100_points(e1_t0, e1_t1, e1_n2_0, e1_n2_1, e1_sf2)[1]
    size_df_1 = pd.DataFrame({'time': time, 'size1': size1, 'size2': size2})
    
    e2_t0 = row['t4']
    e2_t1 = row['t3'] + row['t4']
    e2_n1_0 = row['nu31']
    e2_n1_1 = row['nu21']
    e2_sf1 = pop1_sfs[1]
    time, size1 = get_100_points(e2_t0, e2_t1, e2_n1_0, e2_n1_1, e2_sf1)
    e2_n2_0 = row['nu32']
    e2_n2_1 = row['nu22']
    e2_sf2 = pop2_sfs[1]
    size2 = get_100_points(e2_t0, e2_t1, e2_n2_0, e2_n2_1, e2_sf2)[1]
    size_df_2 = pd.DataFrame({'time': time, 'size1': size1, 'size2': size2})
    
    e3_t0 = row['t3'] + row['t4']
    e3_t1 = row['t2'] + row['t3'] + row['t4']
    e3_n1_0 = row['nu21']
    e3_n1_1 = row['nu11_1']
    e3_sf1 = pop1_sfs[2]
    time, size1 = get_100_points(e3_t0, e3_t1, e3_n1_0, e3_n1_1, e3_sf1)
    e3_n2_0 = row['nu22']
    e3_n2_1 = row['nu11_2']
    e3_sf2 = pop2_sfs[2]
    size2 = get_100_points(e3_t0, e3_t1, e3_n2_0, e3_n2_1, e3_sf2)[1]
    size_df_3 = pd.DataFrame({'time': time, 'size1': size1, 'size2': size2})

    size_df = pd.concat([size_df_1, size_df_2, size_df_3]).reset_index(drop=True)
    size_df['boot'] = boot_num
    return(size_df)

def get_size_df_boot(boot, graph):
    """
    Get the size DataFrame for the whole bootstraps
    """
    size_df_list = []
    for i in range(len(boot)):
        row = boot.iloc[i]
        size_df = get_size_df_row(row, boot_num=i, graph=graph)
        size_df_list.append(size_df)
    return pd.concat(size_df_list).reset_index(drop=True)

log_ticks = [
            100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
            2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
            20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000,
            200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000
]
log_labels1 = [
            '0.1', '', '', '', '', '', '', '', '', '1',
            '', '', '', '', '', '', '', '', '10',
            '', '', '', '', '', '', '', '', '100',
            '', '', '', '', '', '', '', '', '1000'
]

log_labels2 = [
            '', '', '', '', '', '', '', '', '', '1',
            '', '', '', '', '', '', '', '', '10',
            '', '', '', '', '', '', '', '', '100',
            '', '', '', '', '', '', '', '', '1000'
]

fig, axs = plt.subplots(nrows=5, ncols=2, figsize=(5*2, 4.5*2))
fig.subplots_adjust(hspace=0.2, wspace=0.07)

models = ['12_9', '6_2', '20_7', '12_6', '18_7', '18_10']
pairs = ['lpa-wel', 'lpa-wel', 'lpa-wel', 'lpa-sel', 'lpa-sel', 'lpa-sel']


for n in range(6):
    
    # Get the column (first 3 model on the left, last 3 on the right)
    if n in [0, 1, 2]:
        col = 0
    elif n in [3, 4, 5]:
        col = 1
    
    colors = {
        0: '#3f007d', 1: '#41ab5d', 2: '#99000d',
        3: '#3f007d', 4: '#41ab5d', 5: '#99000d',
    }

    # load data
    pop_pair = pairs[n]
    model = models[n]
    model_name = f'{pop_pair}_{model}'
    graph = demes.load(f'data/demographic_inference/{pop_pair}_best_yamls/{model_name}_final_best_model.yaml')
    boot = pd.read_csv(f'data/demographic_inference/{pop_pair}_CI/{model_name}/result_table_converted.csv')
    size_df = get_size_df_boot(boot, graph)
    
    # divergence time
    sns.stripplot(
        x=[n for _ in range(100)],
        y=(boot.t4 + boot.t3 + boot.t2)*5, jitter=0.15,
        ax=axs[0,col], color=colors[n], alpha=0.5
    )
    
    # pop1 size
    for i in range(100):
        df = size_df[size_df['boot'] == i]
        axs[1,col].plot(np.array(df.time)*5, np.array(df.size1), color=colors[n], alpha=0.2)
    
    
    # pop2 size
    for i in range(100):
        df = size_df[size_df['boot'] == i]
        axs[2,col].plot(np.array(df.time)*5, np.array(df.size2), color=colors[n], alpha=0.2)
    

    # migration rate
    for row in boot.iterrows():
        row = row[1]
        x = [0, (row["t4"]), (row["t4"] + row["t3"]), (row["t4"] + row["t3"] + row["t2"])]
        y1 = [row["m4_12"], row["m3_12"], row["m2_12"]]
        y2 = [row["m4_21"], row["m3_21"], row["m2_21"]]
        # axs[3,col].stairs(np.array(y1) * 1e6, np.array(x) * 5, baseline = None, color = colors[n], alpha = 0.2) # recieving pop1 - double check
        # axs[4,col].stairs(np.array(y2) * 1e6, np.array(x) * 5, baseline = None, color = colors[n], alpha = 0.2) # recieving pop2 - double check
        axs[3,col].stairs(np.array(y1), np.array(x) * 5, baseline = None, color = colors[n], alpha = 0.2) # recieving pop1 - double check
        axs[4,col].stairs(np.array(y2), np.array(x) * 5, baseline = None, color = colors[n], alpha = 0.2) # recieving pop2 - double check

### x axis stuff ###
# divergence time plot
axs[0, 0].set_xticks(range(3))
axs[0, 0].set_xticklabels(['model 1', 'model 2', 'model 3'])
axs[0, 0].set_xlabel('ILa - ELw\nDemography\n')
axs[0, 0].xaxis.tick_top()
axs[0, 0].xaxis.set_label_position('top')
axs[0, 1].set_xticks(range(3))
axs[0, 1].set_xticklabels(['model 1', 'model 2', 'model 3'])
axs[0, 1].set_xlabel('ILa - ELs\nDemography\n')
axs[0, 1].xaxis.tick_top()
axs[0, 1].xaxis.set_label_position('top')
# rest of the plots
for i in range(1, 5):
    # log scale
    axs[i, 0].set_xlim(1000, 1_000_000)
    axs[i, 0].set_xscale('log')
    axs[i, 1].set_xlim(1000, 1_000_000)
    axs[i, 1].set_xscale('log')
    if i != 4:
        # axs[i, 0].set_xticks(log_ticks)
        axs[i, 0].set_xticklabels([])
        # axs[i, 1].set_xticks(log_ticks)
        axs[i, 1].set_xticklabels([])
    else:
        axs[i, 0].set_xlabel('Time (years before present)')
        # axs[i, 0].set_xticks(log_ticks)
        # axs[i, 0].set_xticklabels(log_labels1)
        axs[i, 1].set_xlabel('Time (years before present)')
        # axs[i, 1].set_xticks(log_ticks)
        # axs[i, 1].set_xticklabels(log_labels2)

### x axis stuff ###
# divergence time plot
axs[0, 0].set_ylabel('Divergence time\n(thousands of years)\n')
axs[0, 0].set_ylim(0, 1_000_000)
axs[0, 0].set_yticks([0, 200_000, 400_000, 600_000, 800_000, 1_000_000])
axs[0, 0].set_yticklabels(['0', '200', '400', '600', '800', '1000'])
axs[0, 1].set_ylabel('')
axs[0, 1].set_ylim(0, 1_000_000)
axs[0, 1].set_yticks([0, 200_000, 400_000, 600_000, 800_000, 1_000_000])
axs[0, 1].set_yticklabels([])
# size plots
axs[1, 0].set_ylabel('EL population\nsize\n')
axs[2, 0].set_ylabel('IL population\nsize\n')
for i in range(1, 3):
    axs[i, 0].set_ylim(1000, 100_000)
    axs[i, 0].set_yscale('log')
    # axs[i, 0].set_yticks([0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000, 22000, 24000])
    # axs[i, 0].set_yticklabels(['', '2', '', '6', '', '10', '', '14', '', '18', '', '22', ''])
    axs[i, 1].set_ylim(1000, 100_000)
    axs[i, 1].set_yscale('log')
    # axs[i, 1].set_yticks([0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000, 22000, 24000])
    axs[i, 1].set_yticklabels([])
# migration rate plots
axs[3, 0].set_ylabel('IL to EL\nmigration rate\n')
axs[4, 0].set_ylabel('EL to IL\nmigration rate\n')
for i in range(3, 5):
    # axs[i, 0].set_ylim(0, 20)
    axs[i, 0].set_ylim(1e-8, 1e-4)
    axs[i, 0].set_yscale('log')
    axs[i, 0].yaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
    axs[i, 0].yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=10))
    # axs[i, 0].set_yticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7])
    # axs[i, 0].set_yticklabels(['', '', '1', '', '2', '', '3', '', '4', '', '5', '', '6', '', '7'])
    # axs[i, 1].set_ylim(0, 20)
    axs[i, 1].set_ylim(1e-8, 1e-4)
    axs[i, 1].set_yscale('log')
    # axs[i, 1].set_yticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7])
    axs[i, 1].yaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
    axs[i, 1].yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=10))
    axs[i, 1].set_yticklabels([])

plt.savefig('plots/demographic_inference/full_demo_plot.pdf', bbox_inches='tight')
