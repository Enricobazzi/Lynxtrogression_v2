import argparse
from numpy import log
import pandas as pd
import demes
import demesdraw
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser(
    description='Draw the best scoring model using demesdraw and plot the parameter distributions, highlighting the confidence interval limits and the optimal value recovered during the gadma2 run'
    )

parser.add_argument("--params", required=True, help="converted_results_table.csv from convert_ls_on_boot_result.py")
parser.add_argument("--confint", required=True, help="confidence_intervals.tsv from gadma-get_confidence_intervals")
parser.add_argument("--best_yaml", required=True, help="best_model.yaml from gadma2")
parser.add_argument("--output", required=True, help="output directory")

args = parser.parse_args()

run = args.best_yaml.split("/")[-1].split("_")
run = run[1] + "_" + run[2]

## first part of the script: draw the model from the best yaml file
graph = demes.load(args.best_yaml)
# graph.demes[0].epochs[0].start_size = 1
# graph.demes[0].epochs[0].end_size = 1
w = demesdraw.utils.separation_heuristic(graph)

if graph.demes[1].name == "wel_lpa":
    positions = dict(lpa=-w/2, wel=w/2, wel_lpa=0, ancestral=0)
    ax = demesdraw.tubes(graph, log_time=True, positions=positions, labels="xticks-legend", scale_bar=True)
elif graph.demes[1].name == "eel_lpa":
    positions = dict(lpa=-w/2, eel=w/2, eel_lpa=0, ancestral=0)
    ax = demesdraw.tubes(graph, log_time=True, positions=positions, labels="xticks-legend", scale_bar=True)
elif graph.demes[1].name == "sel_lpa":
    positions = dict(lpa=-w/2, sel=w/2, sel_lpa=0, ancestral=0)
    ax = demesdraw.tubes(graph, log_time=True, positions=positions, labels="xticks-legend", scale_bar=True)

plt.savefig(f'{args.output}/{graph.demes[1].name}_{run}_best_model_tubes.png', dpi=300)
plt.close()

ax = demesdraw.size_history(graph, log_time=True, log_size=True)
plt.savefig(f'{args.output}/{graph.demes[1].name}_{run}_best_model_sizes.png', dpi=300)
plt.close()

## second part of the script: plot the parameter distributions
params = pd.read_csv(args.params)
confint = pd.read_csv(args.confint, sep="\t",
                      header=None, names=["parameter", "lower", "upper"])

optim = {}
# parse times
t0 = graph.demes[0].epochs[0].end_time
optim["t1"] = t0 - graph.demes[1].epochs[0].end_time
optim["t2"] = t0 - optim["t1"] - graph.demes[2].epochs[0].end_time
optim["t3"] = t0 - optim["t1"] - optim["t2"] - graph.demes[2].epochs[1].end_time
optim["t4"] = t0 - optim["t1"] - optim["t2"] - optim["t3"] - graph.demes[2].epochs[2].end_time
# parse Ne
optim["Nanc"] = graph.demes[0].epochs[0].end_size
optim["nu11"] = graph.demes[1].epochs[0].end_size
optim["nu11_1"] = graph.demes[2].epochs[0].start_size
optim["nu11_2"] = graph.demes[3].epochs[0].start_size
optim["nu21"] = graph.demes[2].epochs[0].end_size
optim["nu22"] = graph.demes[3].epochs[0].end_size
optim["nu31"] = graph.demes[2].epochs[1].end_size
optim["nu32"] = graph.demes[3].epochs[1].end_size
optim["nu41"] = graph.demes[2].epochs[2].end_size
optim["nu42"] = graph.demes[3].epochs[2].end_size
# parse migration rates
m2_21 = [m.rate for m in graph.migrations if m.source == "wel" and m.dest == "lpa" and round(m.start_time, 3) == round(optim["t4"] + optim["t3"] + optim["t2"], 3)]
if m2_21:
    m2_21 = m2_21[0]
else:
    m2_21 = 0
optim["m2_21"] = m2_21
m2_12 = [m.rate for m in graph.migrations if m.source == "lpa" and m.dest == "wel" and round(m.start_time, 3) == round(optim["t4"] + optim["t3"] + optim["t2"], 3)]
if m2_12:
    m2_12 = m2_12[0]
else:
    m2_12 = 0
optim["m2_12"] = m2_12
m3_21 = [m.rate for m in graph.migrations if m.source == "wel" and m.dest == "lpa" and round(m.start_time, 3) == round(optim["t4"] + optim["t3"], 3)]
if m3_21:
    m3_21 = m3_21[0]
else:
    m3_21 = 0
optim["m3_21"] = m3_21
m3_12 = [m.rate for m in graph.migrations if m.source == "lpa" and m.dest == "wel" and round(m.start_time, 3) == round(optim["t4"] + optim["t3"], 3)]
if m3_12:
    m3_12 = m3_12[0]
else:
    m3_12 = 0
optim["m3_12"] = m3_12
m4_21 = [m.rate for m in graph.migrations if m.source == "wel" and m.dest == "lpa" and round(m.start_time, 3) == round(optim["t4"], 3)]
if m4_21:
    m4_21 = m4_21[0]
else:
    m4_21 = 0
optim["m4_21"] = m4_21
m4_12 = [m.rate for m in graph.migrations if m.source == "lpa" and m.dest == "wel" and round(m.start_time, 3) == round(optim["t4"], 3)]
if m4_12:
    m4_12 = m4_12[0]
else:
    m4_12 = 0
optim["m4_12"] = m4_12

# create the output directory
if not os.path.exists(f'{args.output}/{run}_parameter_distribution'):
    os.makedirs(f'{args.output}/{run}_parameter_distribution')

# plot the parameter distributions
for col in params.columns:
    if col == "bootstrap" or col == "Theta":
        continue
    fig, ax = plt.subplots()
    ax.hist(params[col], bins=30, color="lightblue", edgecolor="black", linewidth=0.5)
    ax.axvline(x=optim[col], color="red", linestyle="--", label="Optimal value")
    ax.axvline(x=confint[confint["parameter"] == col]["lower"].values[0], color="black", linestyle="--", label="lower CI")
    ax.axvline(x=confint[confint["parameter"] == col]["upper"].values[0], color="black", linestyle="--", label="upper CI")
    ax.set_xlabel(col)
    ax.set_ylabel("Frequency")
    ax.legend()
    plt.savefig(f'{args.output}/{run}_parameter_distribution/{col}.png', dpi=300)
    plt.close()
