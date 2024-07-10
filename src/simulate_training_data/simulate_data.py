import pandas as pd
import random
import demes
import re
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Create a tbs file for msmodified from a demes graph and confidence intervals')
    parser.add_argument('--demes_yaml', type=str, help='Path to the demes yaml file', required=True)
    parser.add_argument('--confint', type=str, help='Path to the confidence intervals csv file', required=True)
    parser.add_argument('--path_to_msmodified', type=str, help='Path to the msmodified executable', required=True)
    parser.add_argument('--migration', type=str, help='Migration direction: lly-lpa, lpa-lly, bi, none', required=True)
    parser.add_argument('--odir', type=str, help='Output directory', required=True)
    parser.add_argument('--nreps', type=int, help='Number of replicates', default=5000)
    parser.add_argument('--L', type=int, help='Length of the sequence', default=100_000)
    parser.add_argument('--mu', type=float, help='Mutation rate', default=6e-9)
    parser.add_argument('--rec', type=float, help='Recombination rate', default=1.9e-8)
    return parser.parse_args()

def get_graph_from_yaml(demes_yaml):
    """
    Load the demes graph from a yaml file
    """
    graph = demes.load(demes_yaml)
    return graph

def change_pop1_in_graph_based_on_migration(graph, migration):
    """
    Since msmodified and format.py only work when population 1 is receiving migrants (in forward time),
    we need to change the demes graph accordingly, making the receiving population the ancestral one.

    The default demes graph from gadma has the lynx lynx population as population 1. We therefore have to change the
    demes graph order in case migration is from lynx lynx to lynx pardinus (lly-lpa). In this case lpa from demes[3] in the original
    graph will become demes[0] in graph2, and wel which was demes[2] in the original graph will become demes[1] in graph2.
    On the other hand, if the migration is from lynx pardinus to lynx lynx (lpa-lly), wel will become demes[0] and lpa demes[1] in graph2.

    The ancestral populations epochs are added to whichever becomes demes[0] in the new graph.
    """
    def get_epoch_dict(epoch):
        """
        Get a dictionary version of the epoch from the demes graph - needed by graph.Builder.add_deme(epoch=[dict])
        """
        return {
            "start_size": epoch.start_size,
            "end_size": epoch.end_size,
            "end_time": epoch.end_time,
            "size_function": epoch.size_function
        }
    # new graph
    graph2 = demes.Builder()
    
    if migration == "lly-lpa":
        # demes[3] will be the new demes[0] inheriting the epochs from the ancestral populations plus its own epochs
        deme0_epochs = []
        for epoch in graph.demes[0].epochs:
            deme0_epochs.append(get_epoch_dict(epoch))
        for epoch in graph.demes[1].epochs:
            deme0_epochs.append(get_epoch_dict(epoch))
        for epoch in graph.demes[3].epochs:
            deme0_epochs.append(get_epoch_dict(epoch))
        # add the deme to the new graph
        graph2.add_deme(
            name=graph.demes[3].name,
            epochs=deme0_epochs,
            ancestors=None,
            proportions=None,
            start_time=None,
            description=graph.demes[3].description
        )
        # demes[2] will be the new demes[1] inheriting its own epochs from the original graph
        deme1_epochs = []
        for epoch in graph.demes[2].epochs:
            deme1_epochs.append(get_epoch_dict(epoch))
        # add the deme to the new graph
        graph2.add_deme(
            name=graph.demes[2].name,
            epochs=deme1_epochs,
            ancestors=[graph.demes[3].name],
            proportions=graph.demes[2].proportions,
            start_time=graph.demes[2].start_time,
            description=graph.demes[2].description
        )

    elif migration == "lpa-lly" or migration == "bi" or migration == "none":
        # demes[2] will be the new demes[0] inheriting the epochs from the ancestral populations plus its own epochs
        deme0_epochs = []
        for epoch in graph.demes[0].epochs:
            deme0_epochs.append(get_epoch_dict(epoch))
        for epoch in graph.demes[1].epochs:
            deme0_epochs.append(get_epoch_dict(epoch))
        for epoch in graph.demes[2].epochs:
            deme0_epochs.append(get_epoch_dict(epoch))
        # add the deme to the new graph
        graph2.add_deme(
            name=graph.demes[2].name,
            epochs=deme0_epochs,
            ancestors=None,
            proportions=None,
            start_time=None,
            description=graph.demes[2].description
        )
        # while the demes[3] will be the new demes[1] inheriting its own epochs from the original graph
        deme1_epochs = []
        for epoch in graph.demes[3].epochs:
            deme1_epochs.append(get_epoch_dict(epoch))
        # add the deme to the new graph
        graph2.add_deme(
            name=graph.demes[3].name,
            epochs=deme1_epochs,
            ancestors=[graph.demes[2].name],
            proportions=graph.demes[3].proportions,
            start_time=graph.demes[3].start_time,
            description=graph.demes[3].description
        )
    else:
        raise ValueError("Migration direction not contemplated in the change_pop1_in_graph_based_on_migration script")

    return graph2.resolve()

def get_populations(graph):
    """
    Get the populations names from the demes graph
    """
    pop1 = graph.demes[0].name
    pop2 = graph.demes[1].name
    return pop1, pop2

def get_ms_command(graph, N0, samples, mu, L, rec, migration):
    """
    Get the ms command from the demes graph using demes.to_ms
    """
    theta = 4 * N0 * mu * L 
    rho = rec * 4 * N0 * (L - 1)
    mscommand = demes.to_ms(graph, N0=N0, samples=samples)
    mscommand = f'-t {theta} -r {rho} {L} {mscommand}'
    # add migration to the ms command (not the graph) based on the migration direction
    if migration == "lly-lpa" or migration == "lpa-lly":
        t_split = graph.demes[1].start_time
        mig_time = random.uniform(0, t_split / 4)
        mig_time = mig_time / N0 / 4
        mig_prop = random.uniform(0, 0.5)
        mscommand += f' -es {mig_time} 1 {mig_prop} -ej {mig_time} 3 2'
    elif migration == "bi":
        t_split = graph.demes[1].start_time
        mig_time1 = random.uniform(0, t_split / 4)
        mig_time1 = mig_time1 / N0 / 4
        mig_time2 = random.uniform(0, t_split / 4)
        mig_time2 = mig_time2 / N0 / 4
        mig_prop1 = random.uniform(0, 0.5)
        mig_prop2 = random.uniform(0, 0.5)
        if mig_time1 < mig_time2:
            mscommand += f' -es {mig_time1} 1 {mig_prop1} -ej {mig_time1} 3 2 -es {mig_time2} 2 {mig_prop2} -ej {mig_time2} 4 1'
        else:
            mscommand += f' -es {mig_time2} 2 {mig_prop2} -ej {mig_time2} 3 1 -es {mig_time1} 1 {mig_prop1} -ej {mig_time1} 4 2'
    elif migration == "none":
        mscommand += f' -es 0.0001 1 1.0 -ej 0.0001 3 2'
    else:
        raise ValueError("Migration direction not contemplated in the get_ms_command script")
    return mscommand

def get_ms_string(mscommand, path_to_msmodified, size1, size2, nreps, L):
    """
    Get the ms string with tbs instead of parameters to be used in a bash script
    """
    # get absolute path to msmodified
    path_to_msmodified = os.path.abspath(path_to_msmodified)
    # split the mscommand into the different flags
    pattern = r'(?=\s-[a-zA-Z])'
    commands = [f.strip() for f in re.split(pattern, mscommand)]
    # first part of the ms string
    ms_string = f'{path_to_msmodified} {size1 + size2} {nreps}'
    # add the flags with tbs instead of parameters to the ms string
    for command in commands:
        flag = command.split(' ')[0]
        if flag not in ['-t', '-r', '-I', '-n', '-g', '-en', '-eg', '-ej', '-es']:
            raise ValueError(f'Command from demography not contemplated in the script: {flag}')
        if flag == '-t':
            ms_string += f' {flag} tbs'
        elif flag == '-r':
            ms_string += f' {flag} tbs {L}'
        elif flag == '-I':
            ms_string += f' {flag} {command.split(" ")[1]} {size1} {size2}'
        elif flag in ['-n', '-g']:
            ms_string += f' {flag} {command.split(" ")[1]} tbs'
        elif flag in ['-en', '-eg', '-es']:
            ms_string += f' {flag} tbs {command.split(" ")[2]} tbs'
        elif flag in ['-ej']:
            ms_string += f' {flag} tbs {command.split(" ")[2]} {command.split(" ")[3]}'
    # add the last part of the ms string
    ms_string += f' < mig.tbs | tee mig.msOut\n'
    return ms_string

def get_tbs_list(mscommand):
    pattern = r'(?=\s-[a-zA-Z])'
    commands = [f.strip() for f in re.split(pattern, mscommand)]
    tbs_list = []
    for command in commands:
        flag = command.split(' ')[0]
        if flag not in ['-t', '-r', '-I', '-n', '-g', '-en', '-eg', '-ej', '-es']:
            raise ValueError(f'Command from demography not contemplated in the script: {flag}')
        if flag == '-t':
            tbs_list.append(command.split(' ')[1])
        elif flag == '-r':
            tbs_list.append(command.split(' ')[1])
        elif flag in ['-n', '-g']:
            tbs_list.append(command.split(' ')[2])
        elif flag in ['-en', '-eg', '-es']:
            tbs_list.append(command.split(' ')[1])
            tbs_list.append(command.split(' ')[3])
        elif flag in ['-ej']:
            tbs_list.append(command.split(' ')[1])
    return tbs_list

def modify_graph_from_confint(graph, confint):
    """
    Modify the parameters of the demes graph drawing values from the confidence intervals
    """
    def draw_param_from_confint(param, confint):
        """
        Draw a parameter from the confidence intervals
        """
        return random.uniform(confint.loc[confint["parameter"] == param]["low_ci"].values[0],
                              confint.loc[confint["parameter"] == param]["high_ci"].values[0])
    
    params = {}
    for param in confint["parameter"]:
        params[param] = draw_param_from_confint(param, confint)
    
    ### change parameters in the demes graph:
    # the ancestral population (before split) can be changed without checking population name
    # first epoch in the pre-split population:
    graph.demes[0].epochs[0].end_time = params["t1"] + params["t2"] + params["t3"] + params["t4"]
    graph.demes[0].epochs[0].end_size = params["Nanc"]
    # second epoch in the pre-split population
    if graph.demes[0].epochs[1].size_function == "exponential":
        graph.demes[0].epochs[1].start_size = params["Nanc"]
        graph.demes[0].epochs[1].end_size = params["nu11"]
    else:
        graph.demes[0].epochs[1].start_size = params["nu11"]
        graph.demes[0].epochs[1].end_size = params["nu11"]
    graph.demes[0].epochs[1].end_time = params["t2"] + params["t3"] + params["t4"]
    graph.demes[1].start_time = graph.demes[0].epochs[1].end_time

    # pop1 could be 'lpa' or any of ['wel', 'eel', 'sel] based on pop pair analyzed
    # and who is receiving migrants (pop1 always receives migrants)
    # I can adjust which epoch based on if 'lpa' is the first deme or not
    n = 0 # epoch counter for lpa
    m = 0 # epoch counter for lly = ['wel', 'eel', 'sel]
    if graph.demes[0].name == 'lpa':
        n += 2
    else:
        m += 2
    # define the eurasian lynx population regardless of actual population name
    # based on which of the two demes is not 'lpa'
    if graph.demes[0].name != 'lpa' and graph.demes[0].name in ['wel', 'eel', 'sel']:
        lly = graph.demes[0].name
    elif graph.demes[0].name == 'lpa' and graph.demes[1].name in ['wel', 'eel', 'sel']:
        lly = graph.demes[1].name
    else:
        raise ValueError("Population name from graph not contemplated in the modify_graph_from_confint script")
    
    # first epoch post-split in pop1
    if graph[lly].epochs[m].size_function == "exponential":
        graph[lly].epochs[m].start_size = params["nu11_1"]
        graph[lly].epochs[m].end_size = params["nu21"]
    else:
        graph[lly].epochs[m].start_size = params["nu21"]
        graph[lly].epochs[m].end_size = params["nu21"]
    graph[lly].epochs[m].end_time = params["t3"] + params["t4"]
    # first epoch post-split in pop2
    if graph['lpa'].epochs[n].size_function == "exponential":
        graph['lpa'].epochs[n].start_size = params["nu11_2"]
        graph['lpa'].epochs[n].end_size = params["nu22"]
    else:
        graph['lpa'].epochs[n].start_size = params["nu22"]
        graph['lpa'].epochs[n].end_size = params["nu22"]
    graph['lpa'].epochs[n].end_time = params["t3"] + params["t4"]
    # second epoch post-split in pop1
    if graph[lly].epochs[m+1].size_function == "exponential":
        graph[lly].epochs[m+1].start_size = params["nu21"]
        graph[lly].epochs[m+1].end_size = params["nu31"]
    else:
        graph[lly].epochs[m+1].start_size = params["nu31"]
        graph[lly].epochs[m+1].end_size = params["nu31"]
    graph[lly].epochs[m+1].end_time = params["t4"]
    # second epoch post-split in pop2
    if graph['lpa'].epochs[n+1].size_function == "exponential":
        graph['lpa'].epochs[n+1].start_size = params["nu22"]
        graph['lpa'].epochs[n+1].end_size = params["nu32"]
    else:
        graph['lpa'].epochs[n+1].start_size = params["nu32"]
        graph['lpa'].epochs[n+1].end_size = params["nu32"]
    graph['lpa'].epochs[n+1].end_time = params["t4"]
    # third epoch post-split in pop1
    if graph[lly].epochs[m+2].size_function == "exponential":
        graph[lly].epochs[m+2].start_size = params["nu31"]
        graph[lly].epochs[m+2].end_size = params["nu41"]
    else:
        graph[lly].epochs[m+2].start_size = params["nu41"]
        graph[lly].epochs[m+2].end_size = params["nu41"]
    graph[lly].epochs[m+2].end_time = 0
    # third epoch post-split in pop2
    if graph['lpa'].epochs[n+2].size_function == "exponential":
        graph['lpa'].epochs[n+2].start_size = params["nu32"]
        graph['lpa'].epochs[n+2].end_size = params["nu42"]
    else:
        graph['lpa'].epochs[n+2].start_size = params["nu42"]
        graph['lpa'].epochs[n+2].end_size = params["nu42"]
    graph['lpa'].epochs[n+2].end_time = 0

    return graph

def get_sizes(pop):
    """
    Get the amount of samples to simulate based on the population name
    """
    if pop == 'lpa':
        return 44
    elif pop == 'wel':
        return 40
    elif pop == 'eel':
        return 38
    elif pop == 'sel':
        return 24
    else:
        raise ValueError("Population name not contemplated in the get_sizes script")


def main(demes_yaml, confint, path_to_msmodified, migration, nreps, L, mu, rec, odir):
    """
    Main workflow:
        - Load demes graph from yaml file
        - Add a pulse migration
        - Get ms command with a dummy theta (just for the ms string with tbs instead of parameters)
        - Then for each replicate:
            - Modify the demes graph with parameters drawn from the confidence intervals
            - Get the tbs values for the ms command
    """
    # Load the demes graph from a yaml file
    graph = get_graph_from_yaml(demes_yaml)
    # Change the demes graph based on the migration direction
    graph = change_pop1_in_graph_based_on_migration(graph, migration)
    # Add appropriate migration pulse to the demes graph
    # Load the confidence intervals from the csv file
    confint = pd.read_csv(confint)
    # Get the populations names from the demes graph
    pop1, pop2 = get_populations(graph)
    # Generate a dummy ms command (theta is not important) to get the ms command with tbs instead of parameters 
    mscommand = get_ms_command(graph = graph, N0 = 1000, L=L, samples=[get_sizes(pop1), get_sizes(pop2)], mu=mu, rec=rec, migration=migration)
    ms_string = get_ms_string(mscommand, path_to_msmodified, get_sizes(pop1), get_sizes(pop2), nreps, L)
    # move to odir
    os.chdir(odir)
    # write the tbs file
    with open('mig.tbs', 'w') as f:
        for _ in range(nreps):
            graph = modify_graph_from_confint(graph, confint)
            mscommand = get_ms_command(graph = graph, N0 = graph[pop1].epochs[0].end_size, migration = migration,
                                       L = L, samples = [get_sizes(pop1), get_sizes(pop2)], mu = mu, rec = rec)
            f.write(f"{' '.join(get_tbs_list(mscommand))}\n")
    # write the ms string
    os.system(f'echo "{ms_string}" > ms_string.sh\n')
    os.system(ms_string)
    # gzip useful output (mig.msOut and *.anc)
    os.system(f'gzip mig.msOut')
    os.system(f'gzip *.anc')

if __name__ == '__main__':
    args = parse_args()
    main(args.demes_yaml, args.confint, args.path_to_msmodified, args.migration, args.nreps, args.L, args.mu, args.rec, args.odir)
