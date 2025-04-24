# Genes Analyses

## Get gene translation for Assembly genes

To get the "human" code for all genes in the mLynRuf2.2 assembly I use the information from the Phylome analysis (see Lorena).

My final table will have the following columns for each gene:
- ass_name: name in mLynRuf2.2 assembly
- phy_name: name of LYNRU in the Phylome
- cat_name: name of FELCA in the Phylome
- hum_name: name of gene in human

First two columns come from the `protein_names_0264.txt` file:
```
grep "LYNRU" data/genes_analyses/protein_names_0264.txt |
    awk '{print $2, $1}' | tr ' ' '\t' > data/genes_analyses/gene_names_table.tsv
```

To get the list of phy_name to cat_name I can extract all the Phylome names from the `alns_0264` folder containing all the alignments:
```
for file in $(ls data/genes_analyses/alns_0264); do 
    phy_name=($(grep "LYNRU" data/genes_analyses/alns_0264/${file} | cut -d' ' -f1))
    cat_name=($(grep "FELCA" data/genes_analyses/alns_0264/${file} | cut -d' ' -f1))
    echo -e "${phy_name[@]}\t${cat_name[@]}" | tr ' ' ','
done > data/genes_analyses/LYNRU_to_FELCA_phynames.tsv
```

To get the third column from the `LYNRU_to_FELCA_phynames.tsv` file I just generate:
```
rm -f tmp
phy_names=($(cut -f2 data/genes_analyses/gene_names_table.tsv))
for phy_name in ${phy_names[@]}; do
    cat_names=($(grep ${phy_name} data/genes_analyses/LYNRU_to_FELCA_phynames.tsv | cut -f2 | tr ',' '\n' | sort -u))
    echo -e "${phy_name[@]}\t${cat_names[@]}" | tr ' ' ',' >> tmp
done
paste data/genes_analyses/gene_names_table.tsv <(cut -f2 tmp) > tmp2
mv tmp2 data/genes_analyses/gene_names_table.tsv
rm tmp
```

Fourth column comes from the file `lynx_phylome_ids.txt`:
```
rm -f tmp
cat_names=($(cut -f3 data/genes_analyses/gene_names_table.tsv))
for cat_name in ${cat_names[@]}; do
    hum_names=($(grep -f <(echo ${cat_name} | tr ',' '\n') data/genes_analyses/lynx_phylome_ids.txt | cut -f3))
    echo -e "${cat_name[@]}\t${hum_names[@]}" | tr ' ' ',' >> tmp
done

awk 'NR==FNR {map[$1]=$2; next} {print $0, ($3 in map ? map[$3] : "")}' tmp data/genes_analyses/gene_names_table.tsv |
    tr ' ' '\t' > tmp2

mv tmp2 data/genes_analyses/gene_names_table.tsv
rm tmp
```
*NOTE: the above is super slow and inefficient - should be optimized*

To get the set of genes for each set of introgressed windows:
```
for intro in lpa_to_wel lpa_to_sel wel_to_lpa sel_to_lpa wel_and_sel_to_lpa; do
    echo ${intro}
    grep -f data/introgression_scans/genes/${intro}.geneids.txt data/genes_analyses/gene_names_table.tsv |
        cut -f1,4 > data/genes_analyses/${intro}.hum_genes.txt
done
```