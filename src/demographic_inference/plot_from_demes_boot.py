import demes
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# dummy arguments
graphs = [
    demes.load('data/demographic_inference/lpa-wel_best_yamls/lpa-wel_12_9_final_best_model.yaml'),
    demes.load('data/demographic_inference/lpa-wel_best_yamls/lpa-wel_20_7_final_best_model.yaml'),
    demes.load('data/demographic_inference/lpa-wel_best_yamls/lpa-wel_6_2_final_best_model.yaml'),
]
boots = [
    pd.read_csv('/Users/enrico/Documents/Lynxtrogression_v2/data/demographic_inference/lpa-wel_CI/lpa-wel_12_9/result_table_converted.csv'),
    pd.read_csv('/Users/enrico/Documents/Lynxtrogression_v2/data/demographic_inference/lpa-wel_CI/lpa-wel_20_7/result_table_converted.csv'),
    pd.read_csv('/Users/enrico/Documents/Lynxtrogression_v2/data/demographic_inference/lpa-wel_CI/lpa-wel_6_2/result_table_converted.csv'),
]

log_time = True
num_points = 100

def get_df(graph, idx, num_points=100, log_time=True):
    X = []
    Y = []
    for k, epoch in enumerate(graph.demes[idx].epochs):
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

def modify_graph_from_params(graph, params):
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

for i, graph, boot in zip(range(3), graphs, boots):
    for row in boot.iterrows():
        params = row[1].to_dict()
        graph = modify_graph_from_params(graph, params)
        ancestral = get_df(graph, 1, num_points=num_points, log_time=log_time)
        pop1 = get_df(graph, 2, num_points=num_points, log_time=log_time)
        pop2 = get_df(graph, 3, num_points=num_points, log_time=log_time)
        if i == 0:
            #plt.plot(np.array(ancestral.x * 5), np.array(ancestral.y), label='ancestral', color='royalblue', alpha=0.2)
            plt.plot(np.array(pop1.x * 5), np.array(pop1.y), color='tomato', alpha=0.2)
            plt.plot(np.array(pop2.x * 5), np.array(pop2.y), color='limegreen', alpha=0.2)
        elif i == 1:
            #plt.plot(np.array(ancestral.x * 5), np.array(ancestral.y), color='cornflowerblue', alpha=0.2)
            plt.plot(np.array(pop1.x * 5), np.array(pop1.y), color='orangered', alpha=0.2)
            plt.plot(np.array(pop2.x * 5), np.array(pop2.y), color='forestgreen', alpha=0.2)
        else:
            #plt.plot(np.array(ancestral.x * 5), np.array(ancestral.y), color='blue', alpha=0.2)
            plt.plot(np.array(pop1.x * 5), np.array(pop1.y), color='red', alpha=0.2)
            plt.plot(np.array(pop2.x * 5), np.array(pop2.y), color='green', alpha=0.2)

max_ne = 25000
plt.ylim(0, max_ne)
plt.xlim(100, 2e6)
plt.xscale('log')
plt.xlabel('Years ago')
plt.ylabel('Effective population size')
plt.xticks(ticks=[1e2, 1e3, 1e4, 1e5, 1e6], labels=['100', '1k', '10k', '100k', '1M'])
plt.yticks(ticks=[0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000, 22000, 24000], labels=['0', '2k', '4k', '6k', '8k', '10k', '12k', '14k', '16k', '18k', '20k', '22k', '24k'])
# plt.legend()
plt.show()
