"""
consecutively runs Ising_simulation on a starting structure
that randomly loses 5 edges per simulation run. Two additional 
plots are present in this script
    1) percolation (structure) as a function of degree
    2) synchrony distribition as a function of effective lambda

To run this script, you need to specify the Lambda value in the arguments
i.e. ./more_random_attack.py 1.7

"""

import numpy as np
import networkx as nx
import sys, random, os
import matplotlib.pyplot as plt
import csv
import pandas as pd


from Ising_simulation import IsingSystem, parse_args, writeout_spin_states, preprocess


def prob_in_giant_cluster(structure):
    G = nx.from_numpy_array(structure)
    Gcc = max(nx.connected_components(G), key=len)
    return len(Gcc)/len(G)

def main(args,ver):
    filename = 'input_structure/6025360_20250_2_0_density.txt'
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

    output_dir = 'output/random/{0}/Lam{1}'.format(os.path.basename(filename[:-4]), args.Lambda)
    os.makedirs(output_dir) #make tree of directores

    structure = adjacency
    avg_degrees = []
    P_ones = []
    num_iterations = int(total_edges/num_remove)-80
    synchronies = np.zeros((num_iterations, args.num_sweeps+1))
    for i in range(num_iterations):    #results in around 83 iterations
        for j in range(i*num_remove, (i+1)*num_remove):
            structure[rand_xs[j], rand_ys[j]]=0
            structure[rand_ys[j], rand_xs[j]]=0


        ising_system = IsingSystem(size, args.num_sweeps,
                               args.Lambda, np.copy(structure))
        ising_system.run_simulation(args.verbose)

        writeout_spin_states(output_dir, i, ising_system.output)
        avg_degrees.append(2*(total_edges-(i+1)*num_remove)/size)
        P_ones.append(prob_in_giant_cluster(structure))
        synchronies[i] = ising_system.magnetization_timeseries/size

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

#set to true if wanna check out percolation of structure as edges removed
    if True:
        plt.rcParams.update({'font.size': 14})
        plt.title('$\\Lambda_0 = {0}$, random edge removal'.format(args.Lambda))
        plt.scatter(avg_degrees, P_ones, color='r', label='one simulation')

        plt.xlabel('average degree')
        plt.ylabel('probability in the giant cluster')
        plt.tight_layout()
        plt.show()


#set to true if wanna check out how synchrony changes across targeted attack procedure
    if True:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        nbins = 10
        edges_removed = [0, 10, 20, 30]
        for c, z in zip(['r', 'g', 'b', 'y'], edges_removed):
            ys = synchronies[z] 

            hist, bins = np.histogram(ys, bins=nbins, density=True, range=(-1,1))
            xs = (bins[:-1] + bins[1:])/2

            ax.bar(xs, hist, zs=z, zdir='y', width=1/nbins, color=c, ec=c, alpha=0.8)

        first_lambda = (args.Lambda+1)/(2*size)
        total_edges_removal = total_edges - np.array(edges_removed)*num_remove*2
        probs = total_edges_removal/size**2
        lams = first_lambda*probs * size**2    #multipy by size**2 to get out of normalized united
        print(lams/size**2)

        lams = ["{:.1f}".format(number) for number in lams]
    
        ax.set_yticks(edges_removed) 
        ax.set_yticklabels(lams)    #need to do for only a few decimal places, alos use regular units not /N^2
        ax.set_xlabel('synchrony $s$')
        ax.set_zlabel('$P(s)$')
        ax.set_ylabel('connection strength $\lambda$', labelpad=6)  #labelpad default is 4

        plt.tight_layout()
        plt.show()




if __name__ == "__main__":
    main(parse_args(),'default_name')
