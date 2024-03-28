import os
import pandas as pd
import matplotlib.pyplot as plt

# set the working directory
folder_name = 'data/variant_filtering/missing'

# get a list of the populations
populations = [file.split('.')[0] for file in os.listdir(folder_name) if file.endswith('.nmiss.bed')]

# read in the data
bed_list = []
for pop in populations:
    bed_list.append(pd.read_csv(f'{folder_name}/{pop}.nmiss.bed', sep='\t', header=None))

# get a list of value counts for each bed file
vc_list = []
for bed in bed_list:
    vc = pd.DataFrame(bed[3]).value_counts(sort=False).reset_index()
    vc.columns = ['missing_gts', 'count']
    vc['prop_samples'] = vc['missing_gts'] / max(vc['missing_gts'])
    vc['Freq_prop'] = vc['count'] / len(bed)
    vc['Freq_cumsum'] = vc['Freq_prop'].cumsum()
    vc_list.append(vc)
print(vc_list[0].head())

# plot the cumulative frequency of missing genotypes
plt.figure(figsize=(8, 6))
for val_counts, pop in zip(vc_list, populations):
    plt.plot(val_counts['prop_samples'], val_counts['Freq_cumsum'], alpha=0.5, linewidth=1.5, label=pop)
    plt.scatter(val_counts['prop_samples'], val_counts['Freq_cumsum'], alpha=0.5, s=100, marker='o')
plt.xlim(0, 0.5)
plt.xlabel('Proportion of samples with missing data')
plt.ylabel('Proportion of SNPs included')
plt.legend()
plt.grid(True)
plt.savefig(f'{folder_name}/missing_cumsum.png')

# print a summary per population and write a bed file with the SNPs that have more than 20% missing genotypes
for vc, pop, bed in zip(vc_list, populations, bed_list):
    minmiss = vc[vc["prop_samples"]>0.2]["missing_gts"].iloc[0]
    misssnps = (sum(vc["count"].iloc[minmiss+1:]))
    percmiss = round(vc["Freq_cumsum"].iloc[minmiss]*100, 3)
    print(f'{pop}: {misssnps} SNPs have more than 20% (>{minmiss} samples) missing genotypes ({percmiss}% of total SNPs)')
    filtered_bed = bed[bed[3]>minmiss]
    filtered_bed.to_csv(f'{folder_name}/{pop}.miss_filter.bed', sep='\t', header=False, index=False)


