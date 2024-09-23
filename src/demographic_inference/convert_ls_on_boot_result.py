import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description='Convert the results_table.csv from gadma-run_ls_on_boot_data to physical units'
    )

parser.add_argument("--boot_table", required=True, help="result_table.csv from gadma-run_ls_on_boot_data")
parser.add_argument("--mu", required=True, type=float, help="mutation rate")
parser.add_argument("--L", required=True, type=int, help="length of the sequence")

args = parser.parse_args()

# read the bootstrapping table
df = pd.read_csv(args.boot_table)
# calculate theta0
theta0 = 4 * args.mu * args.L

# change name of first column
df = df.rename(columns={df.columns[0]: "bootstrap"})
# Nanc
df["Nanc"] = df["Theta"] / theta0
# convert parameters in physical units
# iterating over columns
for col in df.columns[1:]:
    if col[0:2] == "nu":
        df[col] = df[col] * df["Nanc"]
    if col[0] == "t":
        df[col] = df[col] * 2 * df["Nanc"]
    if col[0] == "m":
        df[col] = df[col] / (2 * df["Nanc"])

# save converted table
df.to_csv(args.boot_table.replace(".csv", "_converted.csv"), index=False)