"""
phenotype file from UKB.

"""


import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, binned_statistic


def plotout(xs, ys):

    rho, pval = spearmanr(xs, ys)

    means, bins, binnums = binned_statistic(xs, ys, statistic='mean', bins=np.arange(45,85,5))
    std, bins, binnums = binned_statistic(xs, ys, statistic='std', bins=np.arange(45,85,5))
    bins = np.diff(bins)/2.+bins[:-1]

    plt.errorbar(bins, means, yerr = std/np.sqrt(len(xs)), fmt='o-',markersize=10, capsize=6, color='m', label='$\\rho=$ {0:.2f} (<1E-300)\nN = {2}'.format(rho, pval, len(xs)))

    dataset_title = 'UK Biobank'

    plt.title(dataset_title, fontsize=18, fontweight='bold')
    plt.ylabel('white matter volume ($\mathrm{cm}^3$)', fontsize=18)
    plt.xlabel('age in years', fontsize=18)

    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(prop={'size':14})
    
    plt.tight_layout()
    plt.show()


def parse_raw_phenotype_file():
    pheno_file = '/shared/home/rostam/percolating_brain/analyze/phenotypes/ukb/too_big_ukb25909.csv' #not available on GitHUB
                                        #access from UK Biobank

    ages = []
    vols = []

    age_code = '21003-2.0'  #age at dMRI scan
    wm_volume_code = '25007-2.0'   #white matter normalized for head size
#    wm_volume_code = '25008-2.0'   #white matter NOT normalized for head size
#    wm_gm_volume_code = '25009-2.0'   #grey+white normalized for head size


    with open(pheno_file, 'r') as rfile:
        columns = next(rfile).split('","')

        for line in rfile:
            words = line.split('","')
            eid = words[0][1:]  #remove first "

            age = words[columns.index(age_code)]
            wm_volume = words[columns.index(wm_volume_code)]

            if age!='' and wm_volume!='':
                ages.append(float(age))
                vols.append(float(wm_volume)/10**3)
    return ages, vols 



if __name__ == '__main__':
    ages, vols = parse_raw_phenotype_file()

    print('Number of individuals parsed: {0}'.format( len(ages)))

    plotout(ages, vols)
