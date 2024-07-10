import random
import demes
import argparse
import pandas as pd
import random

argparser = argparse.ArgumentParser()
argparser.add_argument("--demes_yaml", required=True, help="Path to the demes graph yaml file with the model to simulate")
argparser.add_argument("--confint", required=True, help="Path to the confidence intervals table")
argparser.add_argument("--size1", required=True, help="Number of samples to simulate from population 1 (eurasian lynx)")
argparser.add_argument("--size2", required=True, help="Number of samples to simulate from population 2 (iberian lynx)")
argparser.add_argument("--L", default=100_000, help="Length of the sequence to simulate")
argparser.add_argument("--mu", default=6e-9, help="Mutation rate")
argparser.add_argument("--rec", default=1.9e-8, help="Recombination rate")
argparser.add_argument("--outfile", help="Output file for the ms commands")
argparser.add_argument("--migration", help="Migration direction", choices=["pop1_to_pop2", "pop2_to_pop1", "both_directions", "no_migration"])
argparser.add_argument("--nreps", default=5000, help="Number of repetitions to simulate")

args = argparser.parse_args()


# load model from yaml file and remove migrations
demes_yaml = args.demes_yaml
graph = demes.load(demes_yaml)
graph.migrations = []

# load confidence intervals from csv file
confint = pd.read_csv(args.confint)

# read arguments for the simulation
size1 = args.size1
size2 = args.size2
L = args.L
mu = args.mu
rec = args.rec

# output file
outfile = args.outfile

# migration direction
migration = args.migration
if migration not in ["pop1_to_pop2", "pop2_to_pop1", "both_directions", "no_migration"]:
    raise ValueError("migration must be one of 'pop1_to_pop2', 'pop2_to_pop1', 'both_directions', or 'no_migration'")

# number of repetitions
nreps = args.nreps

# parse the populations from the demes graph:
pop1 = graph.demes[0].name
pop2 = graph.demes[1].name

# function to draw parameters from the confidence intervals
def draw_param_from_confint(param, confint=confint):
    return random.uniform(confint.loc[confint["parameter"] == param]["low_ci"].values[0],
                          confint.loc[confint["parameter"] == param]["high_ci"].values[0])

for n in range(nreps):
    # create a dictionary with the parameters
    params = {}
    for param in confint["parameter"]:
        params[param] = draw_param_from_confint(param)
    
    ### change parameters in the demes graph:

    # two epochs for the ancestral population (before split - pop1)
    # first epoch in the pre-split population:
    graph[pop1].epochs[0].end_time = params["t1"] + params["t2"] + params["t3"] + params["t4"]
    graph[pop1].epochs[0].end_size = params["Nanc"]
    # second epoch in the pre-split population
    if graph[pop1].epochs[1].size_function == "exponential":
        graph[pop1].epochs[1].start_size = params["Nanc"]
        graph[pop1].epochs[1].end_size = params["nu11"]
    else:
        graph[pop1].epochs[1].start_size = params["nu11"]
        graph[pop1].epochs[1].end_size = params["nu11"]
    graph[pop1].epochs[1].end_time = params["t2"] + params["t3"] + params["t4"]
    graph[pop2].start_time = graph[pop1].epochs[1].end_time
    
    # three epochs in each of the two populations after the split:
    # pop1 epochs 2, 3, 4 and pop2 epochs 0, 1, 2
    # first epoch post-split in pop1
    if graph[pop1].epochs[2].size_function == "exponential":
        graph[pop1].epochs[2].start_size = params["nu11_1"]
        graph[pop1].epochs[2].end_size = params["nu21"]
    else:
        graph[pop1].epochs[2].start_size = params["nu21"]
        graph[pop1].epochs[2].end_size = params["nu21"]
    graph[pop1].epochs[2].end_time = params["t3"] + params["t4"]
    # first epoch post-split in pop2
    if graph[pop2].epochs[0].size_function == "exponential":
        graph[pop2].epochs[0].start_size = params["nu11_2"]
        graph[pop2].epochs[0].end_size = params["nu22"]
    else:
        graph[pop2].epochs[0].start_size = params["nu22"]
        graph[pop2].epochs[0].end_size = params["nu22"]
    graph[pop2].epochs[0].end_time = params["t3"] + params["t4"]
    # second epoch post-split in pop1
    if graph[pop1].epochs[3].size_function == "exponential":
        graph[pop1].epochs[3].start_size = params["nu21"]
        graph[pop1].epochs[3].end_size = params["nu31"]
    else:
        graph[pop1].epochs[3].start_size = params["nu31"]
        graph[pop1].epochs[3].end_size = params["nu31"]
    graph[pop1].epochs[3].end_time = params["t4"]
    # second epoch post-split in pop2
    if graph[pop2].epochs[1].size_function == "exponential":
        graph[pop2].epochs[1].start_size = params["nu22"]
        graph[pop2].epochs[1].end_size = params["nu32"]
    else:
        graph[pop2].epochs[1].start_size = params["nu32"]
        graph[pop2].epochs[1].end_size = params["nu32"]
    graph[pop2].epochs[1].end_time = params["t4"]
    # third epoch post-split in pop1
    if graph[pop1].epochs[4].size_function == "exponential":
        graph[pop1].epochs[4].start_size = params["nu31"]
        graph[pop1].epochs[4].end_size = params["nu41"]
    else:
        graph[pop1].epochs[4].start_size = params["nu41"]
        graph[pop1].epochs[4].end_size = params["nu41"]
    graph[pop1].epochs[4].end_time = 0
    # third epoch post-split in pop2
    if graph[pop2].epochs[2].size_function == "exponential":
        graph[pop2].epochs[2].start_size = params["nu32"]
        graph[pop2].epochs[2].end_size = params["nu42"]
    else:
        graph[pop2].epochs[2].start_size = params["nu42"]
        graph[pop2].epochs[2].end_size = params["nu42"]
    graph[pop2].epochs[2].end_time = 0
    
    ## add pulse migration - migration direction is in forward time
    if migration == "pop1_to_pop2":
        mig_time = random.uniform(0, (params["t4"] + params["t3"] + params["t2"]) / 4)
        mig_prop = random.uniform(0, 0.5)
        if n == 0:
            graph._add_pulse(
                sources=[pop1],
                dest=pop2,
                time=mig_time,
                proportions=[mig_prop]
                )
        else:
            graph.pulses[0].time = mig_time
            graph.pulses[0].proportions = [mig_prop]
    if migration == "pop2_to_pop1":
        mig_time = random.uniform(0, (params["t4"] + params["t3"] + params["t2"]) / 4)
        mig_prop = random.uniform(0, 0.5)
        if n == 0:
            graph._add_pulse(
                sources=[pop2],
                dest=pop1,
                time=mig_time,
                proportions=[mig_prop]
                )
        else:
            graph.pulses[0].time = mig_time
            graph.pulses[0].proportions = [mig_prop]
    if migration == "both_directions":
        mig1_time = random.uniform(0, (params["t4"] + params["t3"] + params["t2"]) / 4)
        mig2_time = random.uniform(0, (params["t4"] + params["t3"] + params["t2"]) / 4)
        mig1_prop = random.uniform(0, 0.5)
        mig2_prop = random.uniform(0, 0.5)
        if n == 0:
            graph._add_pulse(
                sources=[pop1],
                dest=pop2,
                time=mig1_time,
                proportions=[mig1_prop]
                )
            graph._add_pulse(
                sources=[pop2],
                dest=pop1,
                time=mig2_time,
                proportions=[mig2_prop]
                )
        else:
            graph.pulses[0].time = mig1_time
            graph.pulses[0].proportions = [mig1_prop]
            graph.pulses[1].time = mig2_time
            graph.pulses[1].proportions = [mig2_prop]
    if migration == "no_migration":
        pass
    
    # convert the demes graph to ms format
    demes_to_ms = demes.to_ms(graph, N0=params["Nanc"], samples=[size1, size2])
    # calculate theta and rho of the Nanc of this repetition (given mu, rec and L)
    theta = 4 * params["Nanc"] * mu * L 
    rho = rec * 4 * params["Nanc"] * (L - 1)
    # print(f'-t {theta} -r {rho} {L} {demes_to_ms}')
    # write ms command to same file with one line for each replicate (nreps)
    with open(outfile, "a") as f:
        f.write(f'-t {theta} -r {rho} {L} {demes_to_ms}\n')
