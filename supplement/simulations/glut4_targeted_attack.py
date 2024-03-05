"""
consecutively runs Ising_simulation on a starting structure
whose nodes lose edges based on their rank order GLUT4

To run this script, you need to specify the Lambda value in the arguments
i.e. ./glut4_targeted_attack.py 1.7

"""

import numpy as np
import networkx as nx
import sys, os
import matplotlib.pyplot as plt
import pandas as pd

from Ising_simulation import IsingSystem, parse_args, preprocess, writeout_spin_states


def get_average_degree(structure):
    G = nx.from_numpy_array(structure)
    degree = [val for (node, val) in G.degree()]
    return np.mean(degree)

def main(args,ver):

    ta_type = 'decreasing'  #decreasing, increasing

    filename = 'input_structure/6025360_20250_2_0_density.txt'
    pre_adjacency = pd.read_csv(filename, delimiter=" ", header=None).values
    adjacency = preprocess(pre_adjacency)

    adjacency[adjacency>0] = 1    #make adjacency matrix
    size = adjacency.shape[0]


    print("# Initializing Mean-Field ising system of size %d" % (
        size))
    print("# Will run at Lam = %.3f for %d sweeps" % (args.Lambda, args.num_sweeps))



#GLUT4 data mapped to Harvard-Oxford atlas
#obtained from process_glut4/parcellate_map.py
    d_output = {0: 721.045842511998, 1: 3046.1657110041015, 2: 1735.8529609240018, 3: 2716.4412877932805, 4: 777.2103802804043, 5: 1767.7023268573387, 6: 1040.218473521524, 7: 237.16037878155828, 8: 249.3981489102173, 9: 2904.4121281948674, 10: 1864.6278189031368, 11: 2571.401997779131, 12: 662.0312029231615, 13: 1774.9643594150205, 14: 1182.1037845479418, 15: 234.5411436087005, 16: 32625.431783008244, 17: 5765.68268375858, 18: 13151.230638586958, 19: 12356.134211035498, 20: 2630.1558109064326, 21: 3322.1717624881508, 22: 18992.44375855883, 23: 13270.936631555374, 24: 1481.6605808276581, 25: 4565.743035351519, 26: 2568.3322293298716, 27: 6593.014460827609, 28: 4901.430041818143, 29: 1732.2215888981816, 30: 4753.850092762357, 31: 4000.281315962297, 32: 16720.879803965483, 33: 6786.9003496279865, 34: 4210.048850053792, 35: 5705.923319229974, 36: 5943.380883881422, 37: 20606.905737155175, 38: 7860.553767334317, 39: 2673.8291101268387, 40: 2125.1125053727073, 41: 3793.4198861246214, 42: 2479.239337739094, 43: 7752.432606897529, 44: 6142.281878748699, 45: 5179.17143151654, 46: 12971.50773750958, 47: 2391.3520767698633, 48: 7884.500157198087, 49: 2748.139690659568, 50: 1816.2215740629429, 51: 7008.42136891717, 52: 1574.7363231561515, 53: 4405.507646230921, 54: 3471.4125482869276, 55: 3597.1631252059174, 56: 1744.8369030396307, 57: 4524.053379384971, 58: 2686.131314092078, 59: 1694.610368288487, 60: 1372.129921705395, 61: 2177.4934861209213, 62: 570.9236556004078, 63: 8911.390419023533} 


    if ta_type=='increasing':
        sorted_x = sorted(d_output.items(), key=lambda kv: kv[1], reverse=False) 
    elif ta_type=='decreasing':
        sorted_x = sorted(d_output.items(), key=lambda kv: kv[1], reverse=True)


    sub_id = os.path.basename(filename)
    sub_id = sub_id[:sub_id.index('_')]
    output_dir = 'output/targeted/{0}/glut4/{1}/Lam{2}'.format(ta_type, sub_id, args.Lambda)
    os.makedirs(output_dir) #make tree of directores

    structure = adjacency
    avg_degrees = []
    for i, (roi, value) in enumerate(sorted_x):

        structure[roi,:]=0
        structure[:,roi]=0

        ising_system = IsingSystem(size, args.num_sweeps,
                               args.Lambda, np.copy(structure))
        ising_system.run_simulation(args.verbose)
        writeout_spin_states(output_dir, i, ising_system.output)
        avg_degrees.append(get_average_degree(structure))

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
