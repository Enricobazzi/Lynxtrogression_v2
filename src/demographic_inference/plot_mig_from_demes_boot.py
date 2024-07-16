import demes
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# dummy arguments
graphs = [
    demes.load('data/demographic_inference/lpa-wel_best_yamls/lpa-wel_12_9_final_best_model.yaml'),
    #demes.load('data/demographic_inference/lpa-wel_best_yamls/lpa-wel_20_7_final_best_model.yaml'),
    #demes.load('data/demographic_inference/lpa-wel_best_yamls/lpa-wel_6_2_final_best_model.yaml'),
]
boots = [
    pd.read_csv('/Users/enrico/Documents/Lynxtrogression_v2/data/demographic_inference/lpa-wel_CI/lpa-wel_12_9/result_table_converted.csv'),
    #pd.read_csv('/Users/enrico/Documents/Lynxtrogression_v2/data/demographic_inference/lpa-wel_CI/lpa-wel_20_7/result_table_converted.csv'),
    #pd.read_csv('/Users/enrico/Documents/Lynxtrogression_v2/data/demographic_inference/lpa-wel_CI/lpa-wel_6_2/result_table_converted.csv'),
]

def add_migs_to_graph(graph):
    """
    Add migrations to all the epochs in graph if they are not included in the original demes file (migration rate = 0)
    """
    deme0_mig_ends = [m.end_time for m in graph.migrations if m.dest == graph.demes[2].name]
    deme1_mig_ends = [m.end_time for m in graph.migrations if m.dest == graph.demes[3].name]
    for epoch in graph.demes[2].epochs:
        if epoch.end_time not in deme0_mig_ends:
            graph.migrations.append(demes.AsymmetricMigration(
                source=graph.demes[3].name,
                dest=graph.demes[2].name,
                rate=0,
                start_time=epoch.start_time,
                end_time=epoch.end_time
                ))
    for epoch in graph.demes[3].epochs:
        if epoch.end_time not in deme1_mig_ends:
            graph.migrations.append(demes.AsymmetricMigration(
                source=graph.demes[2].name,
                dest=graph.demes[3].name,
                rate=0,
                start_time=epoch.start_time,
                end_time=epoch.end_time
                ))
    return graph

def get_mig_df(graph, idx):
    X = []
    Y = []
    for epoch in graph.demes[idx].epochs:
        start_time = epoch.start_time
        end_time = epoch.end_time
        mig = [m.rate for m in graph.migrations if m.start_time == start_time and m.end_time == end_time][0]
        x = np.array([start_time, end_time])
        y = np.array([mig, mig])
        for i in range(len(x)):
            X.append(x[i])
            Y.append(y[i])
    return pd.DataFrame({'x': X, 'y': Y})

def modify_graph_from_params(graph, params):
    # first epoch
    graph.migrations[0].start_time = params["t4"]
    if graph.migrations[0].dest == 'lpa' and graph.migrations[1].dest == 'wel':
        graph.migrations[0].rate = params["m4_21"]
        graph.migrations[1].rate = params["m4_12"]
    else:
        graph.migrations[0].rate = params["m4_12"]
        graph.migrations[1].rate = params["m4_21"]
    # second epoch
    graph.migrations[2].start_time = params["t3"] + params["t4"]
    graph.migrations[2].end_time = params["t4"]
    if graph.migrations[2].dest == 'lpa' and graph.migrations[3].dest == 'wel':
        graph.migrations[2].rate = params["m3_21"]
        graph.migrations[3].rate = params["m3_12"]
    else:
        graph.migrations[2].rate = params["m3_12"]
        graph.migrations[3].rate = params["m3_21"]
    # third epoch
    graph.migrations[4].start_time = params["t2"] + params["t3"] + params["t4"]
    graph.migrations[4].end_time = params["t3"] + params["t4"]
    if graph.migrations[4].dest == 'lpa' and graph.migrations[5].dest == 'wel':
        graph.migrations[4].rate = params["m2_21"]
        graph.migrations[5].rate = params["m2_12"]
    else:
        graph.migrations[4].rate = params["m2_12"]
        graph.migrations[5].rate = params["m2_21"]
    
    return graph

for i, graph, boot in zip(range(3), graphs, boots):
    graph = add_migs_to_graph(graph)
    for row in boot.iterrows():
        params = row[1].to_dict()
        graph = modify_graph_from_params(graph, params)
        pop1 = get_mig_df(graph, 2)
        pop2 = get_mig_df(graph, 3)
        if i == 0:
            #plt.plot(np.array(ancestral.x * 5), np.array(ancestral.y), label='ancestral', color='royalblue', alpha=0.2)
            plt.plot(np.array(pop1.x * 5), np.array(pop1.y), color='red', alpha=0.2)
            plt.plot(np.array(pop2.x * 5), np.array(pop2.y), color='blue', alpha=0.2)
        elif i == 1:
            #plt.plot(np.array(ancestral.x * 5), np.array(ancestral.y), color='cornflowerblue', alpha=0.2)
            plt.plot(np.array(pop1.x * 5), np.array(pop1.y), color='blue', alpha=0.2)
            plt.plot(np.array(pop2.x * 5), np.array(pop2.y), color='blue', alpha=0.2)
        else:
            #plt.plot(np.array(ancestral.x * 5), np.array(ancestral.y), color='blue', alpha=0.2)
            plt.plot(np.array(pop1.x * 5), np.array(pop1.y), color='green', alpha=0.2)
            plt.plot(np.array(pop2.x * 5), np.array(pop2.y), color='green', alpha=0.2)

plt.xlim(100, 2e6)
plt.xscale('log')
plt.xlabel('Years ago')
plt.yscale('log')
plt.ylabel('Migration rate')
plt.show()
