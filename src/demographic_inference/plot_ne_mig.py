"""
"""
import demes
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

model = '12_6'
pop_pair = 'lpa-sel'
model_name = f'{pop_pair}_{model}'
graph = demes.load(f'data/demographic_inference/{pop_pair}_best_yamls/{model_name}_final_best_model.yaml')
boot = pd.read_csv(f'data/demographic_inference/{pop_pair}_CI/{model_name}/result_table_converted.csv')


log_time = True
num_points = 100

def modify_graph_from_params(graph, params):
    """
    modify a demes graph object with the parameters of a demographic model
    """
    # ancestral population sizes
    graph.demes[0].epochs[0].end_time = params["t1"] + params["t2"] + params["t3"] + params["t4"]
    graph.demes[0].epochs[0].end_size = params["Nanc"]
    # second epoch in the pre-split population
    if graph.demes[1].epochs[0].size_function == "exponential":
        graph.demes[1].epochs[0].start_size = params["Nanc"]
        graph.demes[1].epochs[0].end_size = params["nu11"]
    else:
        graph.demes[1].epochs[0].start_size = params["nu11"]
        graph.demes[1].epochs[0].end_size = params["nu11"]
    graph.demes[1].epochs[0].end_time = params["t2"] + params["t3"] + params["t4"]
    # start time of the split
    graph.demes[2].start_time = graph.demes[1].epochs[0].end_time
    graph.demes[3].start_time = graph.demes[1].epochs[0].end_time
    # first epoch post-split in pop1
    if graph.demes[2].epochs[0].size_function == "exponential":
        graph.demes[2].epochs[0].start_size = params["nu11_1"]
        graph.demes[2].epochs[0].end_size = params["nu21"]
    else:
        graph.demes[2].epochs[0].start_size = params["nu21"]
        graph.demes[2].epochs[0].end_size = params["nu21"]
    graph.demes[2].epochs[0].start_time = params["t2"] + params["t3"] + params["t4"]
    graph.demes[2].epochs[0].end_time = params["t3"] + params["t4"]
    # first epoch post-split in pop2
    if graph.demes[3].epochs[0].size_function == "exponential":
        graph.demes[3].epochs[0].start_size = params["nu11_2"]
        graph.demes[3].epochs[0].end_size = params["nu22"]
    else:
        graph.demes[3].epochs[0].start_size = params["nu22"]
        graph.demes[3].epochs[0].end_size = params["nu22"]
    graph.demes[3].epochs[0].start_time = params["t2"] + params["t3"] + params["t4"]
    graph.demes[3].epochs[0].end_time = params["t3"] + params["t4"]
    # second epoch post-split in pop1
    if graph.demes[2].epochs[1].size_function == "exponential":
        graph.demes[2].epochs[1].start_size = params["nu21"]
        graph.demes[2].epochs[1].end_size = params["nu31"]
    else:
        graph.demes[2].epochs[1].start_size = params["nu31"]
        graph.demes[2].epochs[1].end_size = params["nu31"]
    graph.demes[2].epochs[1].start_time = params["t3"] + params["t4"]
    graph.demes[2].epochs[1].end_time = params["t4"]
    # second epoch post-split in pop2
    if graph.demes[3].epochs[1].size_function == "exponential":
        graph.demes[3].epochs[1].start_size = params["nu22"]
        graph.demes[3].epochs[1].end_size = params["nu32"]
    else:
        graph.demes[3].epochs[1].start_size = params["nu32"]
        graph.demes[3].epochs[1].end_size = params["nu32"]
    graph.demes[3].epochs[1].start_time = params["t3"] + params["t4"]
    graph.demes[3].epochs[1].end_time = params["t4"]
    # third epoch post-split in pop1
    if graph.demes[2].epochs[2].size_function == "exponential":
        graph.demes[2].epochs[2].start_size = params["nu31"]
        graph.demes[2].epochs[2].end_size = params["nu41"]
    else:
        graph.demes[2].epochs[2].start_size = params["nu41"]
        graph.demes[2].epochs[2].end_size = params["nu41"]
    graph.demes[2].epochs[2].start_time = params["t4"]
    graph.demes[2].epochs[2].end_time = 0
    # third epoch post-split in pop2
    if graph.demes[3].epochs[2].size_function == "exponential":
        graph.demes[3].epochs[2].start_size = params["nu32"]
        graph.demes[3].epochs[2].end_size = params["nu42"]
    else:
        graph.demes[3].epochs[2].start_size = params["nu42"]
        graph.demes[3].epochs[2].end_size = params["nu42"]
    graph.demes[3].epochs[2].start_time = params["t4"]
    graph.demes[3].epochs[2].end_time = 0
    
    return graph

def get_df(graph, idx, num_points=100, log_time=True):
    """
    get x and y values for plotting from a particular deme (idx) in a demes graph
    """
    X = []
    Y = []
    for epoch in graph.demes[idx].epochs:
        start_time = epoch.start_time
        if np.isinf(start_time):
            continue
        end_time = epoch.end_time
        if log_time:
            end_time = max(1, end_time)
        if epoch.size_function == "constant":
            x = np.array([start_time, end_time])
            y = np.array([epoch.start_size, epoch.end_size])
        elif epoch.size_function == "exponential":
            x = np.linspace(start_time, end_time, num=num_points)
            dt = np.linspace(0, 1, num=num_points)
            r = np.log(epoch.end_size / epoch.start_size)
            y = epoch.start_size * np.exp(r * dt)
        for i in range(len(x)):
            X.append(x[i])
            Y.append(y[i])
    return pd.DataFrame({'x': X, 'y': Y})

def get_xticks(model_name):
    if model_name == 'lpa-wel_12_9':
        return [1e1, 1e2, 1e3, 1e4, 1e5, 1e6], ['10', '100', '1k', '10k', '100k', '1M']
    else:
        return [1e1, 1e2, 1e3, 1e4, 1e5, 1e6], ['10', '100', '1k', '10k', '100k', '1M']
    
def get_yticks(model_name):
    if model_name == 'lpa-wel_12_9':
        return [0, 2000, 4000, 6000, 8000, 10000], ['0', '2k', '4k', '6k', '8k', '10k']
    else:
        return [0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000, 22000, 24000], ['0', '2k', '4k', '6k', '8k', '10k', '12k', '14k', '16k', '18k', '20k', '22k', '24k']

def plot_ne(graph, boot, opng, num_points=100, log_time=True):
    """
    plot the effective population size over time
    """
    for row in boot.iterrows():
        params = row[1].to_dict()
        graph = modify_graph_from_params(graph, params)
        pop1 = get_df(graph, 2, num_points=num_points, log_time=log_time)
        pop2 = get_df(graph, 3, num_points=num_points, log_time=log_time)
        xticks, xlabels = get_xticks(model_name)
        yticks, ylabels = get_yticks(model_name)
        plt.plot(np.array(pop1.x * 5), np.array(pop1.y), color = 'forestgreen', alpha = 0.7)
        plt.plot(np.array(pop2.x * 5), np.array(pop2.y), color = 'orangered', alpha = 0.7)
    plt.xscale('log')
    plt.xlabel('Years ago')
    plt.ylabel('Effective population size')
    plt.xticks(ticks=xticks, labels=xlabels)
    plt.yticks(ticks=yticks, labels=ylabels)
    plt.title(f'Ne through time in {model_name}')
    plt.ylim(0, 25000)
    plt.gcf().set_size_inches(8, 4)
    plt.savefig(opng, dpi=300)

plot_ne(graph, boot, f'plots/demographic_inference/{model_name}.ne_vs_time.pdf')
plt.clf()

for row in boot.iterrows():
    row = row[1]
    x = [0, (row["t4"]) * 5, (row["t4"] + row["t3"]) * 5, (row["t4"] + row["t3"] + row["t2"]) * 5]
    y1 = [row["m4_12"], row["m3_12"], row["m2_12"]]
    y2 = [row["m4_21"], row["m3_21"], row["m2_21"]]
    plt.stairs(y1, x, baseline = None, color = "forestgreen", alpha = 0.7) # recieving pop1 - double check
    plt.stairs(y2, x, baseline = None, color = "orangered", alpha = 0.7) # recieving pop2 - double check
plt.xscale("log")
# plt.xlim(3, 1e6)
plt.yscale("log")
plt.xlabel('Years ago')
plt.ylabel('Migration rate')
plt.title(f'Migration rate through time in {model_name}')
plt.xticks(ticks=get_xticks(model_name)[0], labels=get_xticks(model_name)[1])
plt.gcf().set_size_inches(8, 4)
plt.savefig(f'plots/demographic_inference/{model_name}.mig_vs_time.pdf', dpi=300)