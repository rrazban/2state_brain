"""
consecutively runs Ising_simulation on a starting structure
that loses edges based on their rank order tract density or
tract length

To run this script, you need to specify the Lambda value in the arguments
i.e. ./targeted_attack.py 1.7

"""

import numpy as np
import sys, os
import matplotlib.pyplot as plt
import pandas as pd

from Ising_simulation import IsingSystem, parse_args, preprocess, writeout_spin_states
from glut4_targeted_attack import get_average_degree


def main(args,ver):
    which = 'density'   #density, length
    ta_type = 'decreasing'  #decreasing, increasing

    filename = 'input_structure/6025360_20250_2_0_{0}.txt'.format(which)
    pre_adjacency = pd.read_csv(filename, delimiter=" ", header=None).values
    weighted_graph = preprocess(pre_adjacency)

#    adjacency[weighted_graph>0] = 1     #don't binarize! cuz studying weighted graph treatment 
    size = weighted_graph.shape[0]


    print("# Initializing Mean-Field ising system of size %d" % (
        size))
    print("# Will run at Lam = %.3f for %d sweeps" % (args.Lambda, args.num_sweeps))


    if ta_type=='increasing':
        thresholds = list(np.logspace(-1, 3, num=100))    #for increasing targeted attack, fully structure not explored but activity is 
    elif ta_type=='decreasing':
        thresholds = list(reversed(np.logspace(-1, 3, num=100)))  #see if can make bounds the same as increasing

    sub_id = os.path.basename(filename)
    sub_id = sub_id[:sub_id.index('_')]
    output_dir = 'output/targeted/{0}/{1}/{2}/Lam{3}'.format(ta_type, which, sub_id, args.Lambda)
    os.makedirs(output_dir) #make tree of directores

    structure = weighted_graph 
    avg_degrees = []

    for i, lim in enumerate(thresholds):       
        if ta_type=='increasing':
            structure[structure<lim] = 0 
        elif ta_type=='decreasing':
            structure[structure>lim] = 0 

        adjacency = np.copy(structure)
        adjacency[adjacency>0] = 1

        ising_system = IsingSystem(size, args.num_sweeps,
                               args.Lambda, adjacency)
        ising_system.run_simulation(args.verbose)
        writeout_spin_states(output_dir, i, ising_system.output)
        avg_degrees.append(get_average_degree(adjacency))

        if False: 
            plt.rcParams.update({'font.size': 14})
            bins = np.arange(-1, 1.1, 0.1)
            plt.hist(ising_system.magnetization_timeseries/size, range=(-1,1))
            plt.title('{0}, $\\Lambda = {1}$'.format(title, args.Lambda))
            plt.xlabel('synchrony')
            plt.ylabel('frequency')
            plt.show()

    with open('{0}/info.csv'.format(output_dir), 'w') as wfile:   #link filename with number of edges left 
        wfile.write('iteration,avg_degree\n')
        for c, cv in enumerate(avg_degrees):
            wfile.write('{0},{1:.2f}\n'.format(c, cv))


if __name__ == "__main__":
    main(parse_args(),'default_name')
