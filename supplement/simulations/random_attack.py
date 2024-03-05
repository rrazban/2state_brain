"""
consecutively runs Ising_simulation on a starting structure
that randomly loses 5 edges per simulation run

To run this script, you need to specify the Lambda value in the arguments
i.e. ./random_attack.py 1.7

"""

import numpy as np
import sys, random, os
import matplotlib.pyplot as plt
import csv
import pandas as pd


from Ising_simulation import IsingSystem, parse_args, writeout_spin_states, preprocess



def main(args,ver):
    filename = 'input_structure/4712851_20250_2_0_density.txt'
    pre_adjacency = pd.read_csv(filename, delimiter=" ", header=None).values
    adjacency = preprocess(pre_adjacency)

    adjacency[adjacency>0] = 1    #make adjacency matrix
    size = adjacency.shape[0]
	

    print("# Initializing Mean-Field ising system of size %d" % (
        size))
    print("# Will run at Lam = %.3f for %d sweeps" % (args.Lambda, args.num_sweeps))


   
    M_tri = np.copy(adjacency)
    M_tri[np.triu_indices(size)] = 0   #find unique edges   
    xs, ys = np.nonzero(M_tri)    #will remove both forward and reverse edges at once
    total_edges = len(xs)
    indis = list(range(total_edges))
    random.shuffle(indis)   #randomize to choose which edges to remove first

    rand_xs = xs[indis]
    rand_ys = ys[indis]
    
    num_remove = 5  #technically *2 cuz remove both front and back connections

    sub_id = os.path.basename(filename)
    sub_id = sub_id[:sub_id.index('_')]
    output_dir = 'output/random/{0}/Lam{1}'.format(sub_id, args.Lambda)
    os.makedirs(output_dir) #make tree of directores

    structure = adjacency
    avg_degrees = []
    for i in range(int(total_edges/num_remove)-80):    #results in around 83 iterations
        for j in range(i*num_remove, (i+1)*num_remove):
            structure[rand_xs[j], rand_ys[j]]=0
            structure[rand_ys[j], rand_xs[j]]=0


        ising_system = IsingSystem(size, args.num_sweeps,
                               args.Lambda, np.copy(structure))
        ising_system.run_simulation(args.verbose)

        writeout_spin_states(output_dir, i, ising_system.output)
        avg_degrees.append(2*(total_edges-(i+1)*num_remove)/size)

#set to true to check out simulated synchrony distribution 
#needs to be at least roughly symmetric to satistify equilibrium
        if False: 
            plt.rcParams.update({'font.size': 14})
            bins = np.arange(-1, 1.1, 0.1)
            plt.hist(ising_system.magnetization_timeseries/size, range=(-1,1))
            plt.title('$\\Lambda_0 = {0}$, {1} edges removed'.format(args.Lambda, (i+1)*num_remove))
            plt.xlabel('synchrony')
            plt.ylabel('frequency')
            plt.show()
 

    with open('{0}/info.csv'.format(output_dir), 'w') as wfile:   #link filename with number of edges left 
        wfile.write('iteration,avg_degree\n')
        for c, cv in enumerate(avg_degrees):
            wfile.write('{0},{1:.2f}\n'.format(c, cv))


if __name__ == "__main__":
    main(parse_args(),'default_name')
