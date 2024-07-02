## plot the parameter distributions and write the median and 95% CI to a table
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse

# parse the arguments
parser = argparse.ArgumentParser(
    "Plot the parameter distributions and write the median and 95% CI to a table"
)
parser.add_argument("--csv", required=True, help="Path to the csv file with the bootstrapped parameters")
parser.add_argument("--model", required=True, help="Model name")
parser.add_argument("--pop_pair", required=True, help="Population pair")
args = parser.parse_args()

csv = args.csv
model = args.model
pop_pair = args.pop_pair

# load the data
boot_table = pd.read_csv(csv)

# plot distributions
if not os.path.exists(f'plots/demographic_inference/{pop_pair}_{model}_parameter_distribution'):
    os.makedirs(f'plots/demographic_inference/{pop_pair}_{model}_parameter_distribution')

for col in boot_table.columns:
    if col == "bootstrap" or col == "Theta":
        continue
    else:
        ## calculate the median and 95% CI
        # median = median of the column
        median = boot_table[col].median()
        # low CI = 2.5th percentile
        low_ci = boot_table[col].quantile(0.025)
        # high CI = 97.5th percentile
        high_ci = boot_table[col].quantile(0.975)
        # plot distribution + median + 95% CI
        fig, ax = plt.subplots()
        ax.hist(boot_table[col], bins=30, color='lightblue')
        ax.axvline(median, color='red', linestyle='dashed', linewidth=0.75, label='Median')
        ax.axvline(low_ci, color='black', linestyle='dashed', linewidth=0.75, label='CI 2.5%')
        ax.axvline(high_ci, color='black', linestyle='dashed', linewidth=0.75, label='CI 97.5%')
        ax.set_xlabel(col)
        ax.set_ylabel("Frequency")
        ax.legend()
        plt.savefig(f'plots/demographic_inference/{pop_pair}_{model}_parameter_distribution/{col}.median.png', dpi=300)
        plt.close()

# create a table with the median and 95% CI
ci_table = pd.DataFrame(columns=["parameter", "median", "low_ci", "high_ci"])
params = boot_table.columns[1:].values
median = [boot_table[col].median() for col in params]
low_ci = [boot_table[col].quantile(0.025) for col in params]
high_ci = [boot_table[col].quantile(0.975) for col in params]
ci_table["parameter"] = params
ci_table["median"] = median
ci_table["low_ci"] = low_ci
ci_table["high_ci"] = high_ci
ci_table.to_csv(f'data/demographic_inference/{pop_pair}_CI/{pop_pair}.{model}.CI.csv', index=False)
