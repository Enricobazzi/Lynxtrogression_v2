import random
import demes
import argparse
import pandas as pd
import random
import os

# parser = argparse.ArgumentParser("Create ms command for a given demes graph drawing parameters from a confidence intervals table")
# 
# parser.add_argument("demes_yaml", help="Path to the demes graph yaml file")
# parser.add_argument("confidence_intervals", help="Path to the confidence intervals table")
# parser.add_argument("size1", help="Number of samples to simulate from population 1")
# parser.add_argument("size2", help="Number of samples to simulate from population 2")
# 
# args = parser.parse_args()
# 
# graph = demes.load(args.demes_yaml)
# confint = pd.read_csv(args.confidence_intervals, sep="\t",
#                       header=None, names=["parameter", "lower", "upper"])

graph = demes.load("/Users/enrico/Documents/Lynxtrogression_v2/data/demographic_inference/lpa-wel_best_yamls/lpa-wel_20_7_final_best_model.yaml")
confint = pd.read_csv("/Users/enrico/Documents/Lynxtrogression_v2/data/demographic_inference/lpa-wel_CI/lpa-wel_20_7/confidence_intervals.tsv", sep="\t",
                      header=None, names=["parameter", "lower", "upper"])
size1 = 40
size2 = 44
L = 100_000
mu = 6e-9
rec = 1.9e-8 

# parse the populations from the demes graph:
# the first two (graph.demes[0] and graph.demes[1]) are the ancestral populations (before split)
pop1 = graph.demes[2].name
pop2 = graph.demes[3].name

# remove migrations from the demes graph
graph.migrations = []

# function to draw parameters from the confidence intervals
def draw_param_from_confint(param, confint=confint):
    return random.uniform(confint.loc[confint["parameter"] == param]["lower"].values[0],
                          confint.loc[confint["parameter"] == param]["upper"].values[0])

# create a dictionary with the parameters
params = {}
for param in confint["parameter"]:
    params[param] = draw_param_from_confint(param)

# change parameters in the demes graph:

# ancestral end time
graph["ancestral"].epochs[0].end_time = params["t1"] + params["t2"] + params["t3"] + params["t4"]
# ancestral end size
graph["ancestral"].epochs[0].end_size = params["Nanc"]

# pre-split population size
if graph.demes[1].epochs[0].size_function == "exponential":
    graph.demes[1].epochs[0].start_size = params["Nanc"]
    graph.demes[1].epochs[0].end_size = params["nu11"]
else:
    graph.demes[1].epochs[0].start_size = params["nu11"]
    graph.demes[1].epochs[0].end_size = params["nu11"]
# pre-split start time
graph.demes[1].start_time = graph["ancestral"].epochs[0].end_time
# pre-split end time
graph.demes[1].epochs[0].end_time = params["t2"] + params["t3"] + params["t4"]

# population 1 and 2 start time
graph[pop1].start_time = graph.demes[1].epochs[0].end_time
graph[pop2].start_time = graph.demes[1].epochs[0].end_time

# population 1 size epoch 0
if graph[pop1].epochs[0].size_function == "exponential":
    graph[pop1].epochs[0].start_size = params["nu11_1"]
    graph[pop1].epochs[0].end_size = params["nu21"]
else:
    graph[pop1].epochs[0].start_size = params["nu21"]
    graph[pop1].epochs[0].end_size = params["nu21"]
# population 1 end time
graph[pop1].epochs[0].end_time = params["t3"] + params["t4"]

# population 2 size epoch 0
if graph[pop2].epochs[0].size_function == "exponential":
    graph[pop2].epochs[0].start_size = params["nu11_2"]
    graph[pop2].epochs[0].end_size = params["nu22"]
else:
    graph[pop2].epochs[0].start_size = params["nu22"]
    graph[pop2].epochs[0].end_size = params["nu22"]
# population 2 end time
graph[pop2].epochs[0].end_time = params["t3"] + params["t4"]

# population 1 size epoch 1
if graph[pop1].epochs[1].size_function == "exponential":
    graph[pop1].epochs[1].start_size = params["nu21"]
    graph[pop1].epochs[1].end_size = params["nu31"]
else:
    graph[pop1].epochs[1].start_size = params["nu31"]
    graph[pop1].epochs[1].end_size = params["nu31"]
# population 1 end time
graph[pop1].epochs[1].end_time = params["t4"]

# population 2 size epoch 1
if graph[pop2].epochs[1].size_function == "exponential":
    graph[pop2].epochs[1].start_size = params["nu22"]
    graph[pop2].epochs[1].end_size = params["nu32"]
else:
    graph[pop2].epochs[1].start_size = params["nu32"]
    graph[pop2].epochs[1].end_size = params["nu32"]
# population 2 end time
graph[pop2].epochs[1].end_time = params["t4"]

# population 1 size epoch 2
if graph[pop1].epochs[2].size_function == "exponential":
    graph[pop1].epochs[2].start_size = params["nu31"]
    graph[pop1].epochs[2].end_size = params["nu41"]
else:
    graph[pop1].epochs[2].start_size = params["nu41"]
    graph[pop1].epochs[2].end_size = params["nu41"]
# population 1 end time
graph[pop1].epochs[2].end_time = 0

# population 2 size epoch 2
if graph[pop2].epochs[2].size_function == "exponential":
    graph[pop2].epochs[2].start_size = params["nu32"]
    graph[pop2].epochs[2].end_size = params["nu42"]
else:
    graph[pop2].epochs[2].start_size = params["nu42"]
    graph[pop2].epochs[2].end_size = params["nu42"]
# population 2 end time
graph[pop2].epochs[2].end_time = 0

# add pulse migration
mig_time = random.uniform(0, (params["t4"] + params["t3"] + params["t2"]) / 4)
mig_prop = random.uniform(0, 1)

graph._add_pulse(
    sources=[pop1],
    dest=pop2,
    time=mig_time,
    proportions=[mig_prop]
    )

demes_to_ms = demes.to_ms(graph, N0 = params["Nanc"],
                          samples = [0, 0, size1, size2])
nsam = size1 + size2
theta = 4 * params["Nanc"] * mu * L 
rho = rec * 4 * params["Nanc"] * (L - 1)
print(f'msmodified/ms {nsam} 1 -t {theta} -r {rho} {L} {demes_to_ms}')
# os.system(f'msmodified/ms {nsam} 1 -t {theta} -r {rho} {L} {demes_to_ms}')

## TODO: need to remove populations 1 and 2 (ancestrals) and make populations 3 and 4 as 1 and 2.
# Try from the demes graph first - DO I THOUGH? EH? DO I?