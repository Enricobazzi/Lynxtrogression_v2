import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    'Get the confidence intervals for the parameters from the converted optimization table and plot the parameter distributions'
    )

parser.add_argument("--table", required=True, help="converted optimization table")
parser.add_argument("--out", required=True, help="output prefix for the table and plots")
args = parser.parse_args()

# read the converted optimization table
df = pd.read_csv(args.table)

# get the parameters
params = df.columns[4:].to_list()

# empty df to store the confidence intervals
ci = pd.DataFrame(columns=["parameter", "lower", "upper"])

# for each parameter, get the confidence intervals and plot the distribution 
for param in params:
    # get the mean, lower and upper confidence intervals
    lower = df[param].quantile(0.025)
    upper = df[param].quantile(0.975)
    # append to the ci df
    ci = pd.concat([ci, pd.DataFrame({"parameter": [param], "lower": [lower], "upper": [upper]})], ignore_index=True)
    # plot the distribution with the confidence intervals and the real data value
    real_data = df[df["data"] == "real_data"][param].values[0]
    # print(f'{param}: {real_data} [{lower} - {upper}]')
    fig, ax = plt.subplots()
    ax.hist(df[param], bins=30, color="lightblue", edgecolor="black", linewidth=0.5)
    ax.axvline(x=real_data, color="red", linestyle="--", label="real data value")
    ax.axvline(x=lower, color="black", linestyle="--", label="lower CI")
    ax.axvline(x=upper, color="black", linestyle="--", label="upper CI")
    ax.set_xlabel(param)
    ax.set_ylabel("Frequency")
    ax.legend()
    plt.savefig(f'plots/demographic_inference/{args.out}.{param}.png', dpi=300)
    plt.close()

# save the confidence intervals
ci.to_csv(f'data/demographic_inference/{args.out}.ci.csv', index=False)