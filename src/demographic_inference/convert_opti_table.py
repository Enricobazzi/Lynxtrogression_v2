import argparse
import pandas as pd
import dadi
import gadma2_best_dadi_models as models

parser = argparse.ArgumentParser(
    description='Convert the <pop_pair>_<model>.opti_table.csv from the dadi optimizations to physical units'
    )

parser.add_argument("--pop_pair", required=True, help="population pair")
parser.add_argument("--model", required=True, help="model used for the optimization")
parser.add_argument("--mu", required=True, type=float, help="mutation rate")
parser.add_argument("--L", required=True, type=int, help="length of the sequence")

args = parser.parse_args()

# read the bootstrapping table
opti_table = f'data/demographic_inference/{args.pop_pair}_CI/{args.pop_pair}_{args.model}.opti_table.csv'
df = pd.read_csv(opti_table)

# read the data folder
data_folder = f'data/demographic_inference/{args.pop_pair}_bootstrap'

# pts
pts = [50, 60, 70]

# get the model
model_func = dadi.Numerics.make_extrap_log_func(getattr(models, f"model_func_{args.model}"))
# model = model_func(p0p, data.sample_sizes, pts)

# TODO: THIS SHOULD GO IN THE DADI_OPTIMIZATION.PY SCRIPT SO THAT THE TABLE ALREADY HAS THE THETA VALUES
# calculate optimal theta for each row
thetas = []
for data in df["data"].unique():
    if "real" in data:
        name = data
    else:
        name = data.split("_")[1]
    fs = dadi.Spectrum.from_file(f'{data_folder}/{args.pop_pair}_{name}.fs')
    params = df[df["data"] == data].iloc[0, 4:].values
    model = model_func(params, fs.sample_sizes, pts)
    print(f'calculating theta for {data}...', end='\r')
    theta = dadi.Inference.optimal_sfs_scaling(model, fs)
    thetas.append(float(theta))

df["Theta"] = thetas

# calculate theta0
theta0 = 4 * args.mu * args.L

# Nanc
df["Nanc"] = df["Theta"]/theta0
# convert parameters in physical units
# iterating over columns
for col in df.columns[1:]:
    if col[0:2] == "nu":
        df[col] = df[col] * df["Nanc"]
    if col[0] == "t":
        df[col] = df[col] * 2 * df["Nanc"]
    if col[0] == "m" and col != "model":
        df[col] = df[col] / (2 * df["Nanc"])

# save converted table
df.to_csv(opti_table.replace(".csv", ".converted.csv"), index=False)