"""
Code to simulate the mean field Ising model.

To run this script, you need to specify the Lambda value in the arguments
i.e. ./Ising_simulation.py 1.7

"""

import numpy as np
import sys, random, os
import matplotlib.pyplot as plt
import csv
import pandas as pd
import networkx as nx


def identify_critical_point(adjacency):
	#this provides an upper bound but still need to search
	G = nx.from_numpy_array(adjacency)
	degree = [val for (node, val) in G.degree()]
	p_edge = np.mean(degree)/adjacency.shape[0]

	return 1/p_edge -1


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
#    parser.add_argument('size', type=int, metavar='N',
 #                       help="size of system") #no need because specified by diffusion MRI structure
    parser.add_argument('Lambda', type=float, metavar='LAMBDA', help="units of coupling/kT. rescaled such that for a fully connected graph, the critical point is at Lam=0.")
    parser.add_argument('--num-sweeps', type=int, metavar='NUM', default=2500)	#corresponds to the number of measurements
    parser.add_argument('--verbose', action='store_true')
    return parser.parse_args()


def boltzmann_weight(energy, beta):
    return np.exp(energy * beta)

def preprocess(r):
    print_extra_info = False 

    np.fill_diagonal(r, 0)	#if dont remove diagonals, shifts avg degree to the right!

    r2 = r[~np.all(r == 0, axis=1)]
    adjacency_matrix = r2[:, ~np.all(r2 == 0, axis=0)]

    if print_extra_info:
        print(list(np.where(~r.any(axis=1))[0]))    #indices of excluded regions
        print(r.shape)
        print(r2.shape)
        print(adjacency_matrix.shape)

    return adjacency_matrix


class IsingSystem:
    """ 
    the Mean-Field Ising System
    """
    def __init__(self, size, sweeps, Lambda, adjacency):

        self.system = 2*np.random.randint(2, size=size) - 1 #initial spin state randomly set

        self.lambdaa = (Lambda+1)/(2*size) #note that in the traditional Ising model formalism, lambda refers to coupling/temperature ($\lambda$=J/T)
#       self.lambdaa = np.random.normal(loc=self.lambdaa, scale=3*self.beta, size=((size,size)))  #uncomment if want to explore distribution of lambdas 


        self.adjacency = adjacency
        self.times = sweeps+1
        self.size = size

        self.magnetization_timeseries = np.zeros(sweeps+1)
        self.energy_timeseries = np.zeros(sweeps+1)
        self.success_timeseries = np.zeros(sweeps+1, dtype=int)        # track number of accepted moves per sweep:
        self.output = np.zeros((sweeps+1, size))


        """
        The following two variables set the number of steps per sweep and the number of attempted spin flips per step 
        They are important to get the simulation to quickly reach equilibrium
            -we need to reach equilibrium if we want to compare to mathematical equations because they are derived at equilibrium
        Their exact values do not matter, as long as equilibrium is reached
        """
        self.steps_per_sweep = size 
        self.num_attempts = int(0.15*size)#10	#num of spins attempted to be flipped at a given time step

        return

    def calculate_magnetization(self):	#otherwise known as synchrony when taking the mean
       # return np.mean(self.system)
        return np.sum(self.system)

    def calculate_energy(self, spins):
        """
        The full energy function involves calculation of all pairwise
        energy interactions. Assume spins are in contact according to the adjacency matrix 
        """

        energy = self.lambdaa * np.dot(spins, np.dot(self.adjacency, spins))
        return energy       #even though there is double counting of edges, do not divide by 2 to be consistent with fully connected ising model formalism
    
    def calculate_deltaE(self, positions):
        """
        Position should be an index corresponding to a region.
        """
        new_spins = np.copy(self.system)
        new_spins[positions] *= -1

        deltaE = self.calculate_energy(new_spins) - self.calculate_energy(self.system)
        return deltaE
            

    def run_simulation(self, verbose=False):
        """
        Implements a Metropolis-Hastings algorithm to determine whether 
        attempted sets of spin flips are accepted or not
        """


        """initialize output variables"""
        self.output[0] = self.system
        self.magnetization_timeseries[0] = self.calculate_magnetization()
        self.energy_timeseries[0] = self.calculate_energy(self.system)
        self.success_timeseries[0] = 0
        if verbose:
            print("{0:>10} {1:>15} {2:>15} {3:>15}".format(
                "Time", "Magnetization", "Energy", "NumSuccesses"))
            print("{0:10d} {1:15.3f} {2:15.3f} {3:15d}".format(
                0, self.magnetization_timeseries[0],
                self.energy_timeseries[0], self.success_timeseries[0]))

        """let the dynamics begin"""
        for time in range(1,self.times):
            num_successes = 0
            for j in range(self.steps_per_sweep):
                selected_spins = np.random.randint(self.size, size=self.num_attempts)

                deltaE = self.calculate_deltaE(selected_spins)
                """Metropolis-Hastings algorithm"""
                if deltaE > 0:
                    self.system[selected_spins] *= -1
                    num_successes += 1
                else:
                    if np.random.rand() < boltzmann_weight(deltaE, 1):
                        self.system[selected_spins] *= -1
                        num_successes += 1

            """record output values"""
            self.magnetization_timeseries[time] = self.calculate_magnetization()
            self.energy_timeseries[time] = self.calculate_energy(self.system)
            self.success_timeseries[time] = num_successes
            self.output[time] = self.system
            if verbose:
                if num_successes>0: 
                	print("{0:10d} {1:15.3f} {2:15.3f} {3:15d}".format(
                    	time, self.magnetization_timeseries[time],
                    	self.energy_timeseries[time], self.success_timeseries[time]))


        return

def writeout_spin_states(dirname, sim_id, output):
	filename = 'sim-{0}.csv'.format(sim_id)
	print('output filename: {0}/{1}'.format(dirname, filename))
	with open("{0}/{1}".format(dirname, filename), "w", newline='') as f:
		writer = csv.writer(f)
		writer.writerows(output) 


def main(args,ver):

    filename = 'input_structure/6025360_20250_2_0_density.txt'   #this structure has 64 regions with at least one connection
    pre_adjacency = pd.read_csv(filename, delimiter=" ", header=None).values
    adjacency = preprocess(pre_adjacency)   #remove those regions with no connections which gets you to 64 connections

    adjacency[adjacency>0] = 1    #make adjacency matrix
    size = adjacency.shape[0]
	
    print("# Initializing Mean-Field ising system with %d nodes" % (
        size))
    print("# Will run at Lam = %.3f for %d sweeps" % (args.Lambda, args.num_sweeps))
#    print('# The critical point of this system is below Lam = %.3f' % (identify_critical_point(adjacency)))


    ising_system = IsingSystem(size, args.num_sweeps,
                               args.Lambda, np.copy(adjacency))
    ising_system.run_simulation(args.verbose)

    sim_id = 0
    output_dir = 'output/individual/{0}/Lam{1}'.format(os.path.basename(filename[:-4]), args.Lambda)
    os.makedirs(output_dir) #make tree of directores
    writeout_spin_states(output_dir, sim_id, ising_system.output)

    if False:   #set to true to check out simulated synchrony distribution 
            plt.rcParams.update({'font.size': 14})
            plt.hist(ising_system.magnetization_timeseries/size, range=(-1,1))
            plt.title('$\\Lambda = {0}$'.format(args.Lambda))
            plt.xlabel('synchrony')
            plt.ylabel('frequency')
            plt.tight_layout()
            plt.show()
    

if __name__ == "__main__":
    main(parse_args(),'default_name')
