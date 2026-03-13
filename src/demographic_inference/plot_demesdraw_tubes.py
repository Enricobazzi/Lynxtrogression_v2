import demes
import demesdraw
import matplotlib.pyplot as plt

pop_pairs = [
    'lpa-wel',
    'lpa-sel',
]

#for pair in pop_pairs:
for pair in pop_pairs:
    if pair == "lpa-wel":
        models = ["6_2", "12_9", "20_7"]
    elif pair == "lpa-sel":
        models = ["12_6", "18_7", "18_10"]
    for model in models:
        graph = demes.load(f'data/demographic_inference/{pair}_best_yamls/{pair}_{model}_final_best_model.yaml')
        # remove ancestral populations for clearer plots
        graph.demes = [d for d in graph.demes if d.name in [pair.split("-")[0], pair.split("-")[1]]]
        graph.demes[0].ancestors = []
        graph.demes[0].proportions = []
        graph.demes[1].ancestors = []
        graph.demes[1].proportions = []
        w = demesdraw.utils.separation_heuristic(graph)
        # set positions
        if pair == "lpa-wel":
            positions = dict(lpa=-w/2, wel=w/2, wel_lpa=0)
        elif pair == "lpa-sel":
            positions = dict(lpa=-w/2, sel=w/2, sel_lpa=0)
        # create and save plot
        fig, ax = plt.subplots(figsize=(6, 4))
        demesdraw.tubes(graph, log_time=True, positions=positions, labels="xticks-legend", scale_bar=False)
        fig = plt.gcf()
        fig.set_size_inches(5, 5)
        fig.savefig(
            f'plots/revision1/tubes/{pair}_{model}.pdf',
            dpi=300,
            format='pdf',
            bbox_inches='tight'
        )
        plt.close(fig)