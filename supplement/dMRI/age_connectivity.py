"""
Plot dMRI properties as a function of ages. Data 
is binned for ease in visualization, however, 
correlations are calculated by considering data 
points

"""


import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import spearmanr, binned_statistic


def plotout(xs, ys, y2s, which):

    if which!='avg_degree':
        ys = np.array(ys)/10**6
        y2s = np.array(y2s)/10**6

    rho, pval = spearmanr(xs, ys)
    rho2, pval2 = spearmanr(xs, y2s)

    bins = np.arange(45, 85, 5)
    means, bins, binnums = binned_statistic(xs, ys, statistic='mean', bins=bins)
    std, bins, binnums = binned_statistic(xs, ys, statistic='std', bins=bins)
    means2, bins2, binnums = binned_statistic(xs, y2s, statistic='mean', bins=bins)
    std2, bins2, binnums = binned_statistic(xs, y2s, statistic='std', bins=bins)
    bins = np.diff(bins)/2.+bins[:-1]

    plt.errorbar(bins, means, yerr = std/np.sqrt(len(xs)), fmt='o-',markersize=10, capsize=6, color='r', label='constant solid angle (Q-Ball)\n$\\rho=$ {0:.2f} ({1:.2E})'.format(rho, pval, len(xs)))
    plt.errorbar(bins, means2, yerr = std2/np.sqrt(len(xs)), fmt='o-',markersize=10, capsize=6, color='b', label='diffusion tensor imaging\n$\\rho=$ {0:.2f} ({1:.2E})'.format(rho2, pval2, len(xs)))

    plt.title('UK Biobank', fontsize=18, fontweight='bold')


    if which=='density':
        #plt.ylabel('total tract {0} (per million counts)'.format(which), fontsize=18)  #overflow
        plt.ylabel('total tract {0} (per million)'.format(which), fontsize=18)
    elif which=='length':
        plt.ylabel('total tract {0} (km)'.format(which), fontsize=18)
    else:
        plt.ylabel('average degree', fontsize=18)


    plt.xlabel('age in years', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    legend = plt.legend(prop={'size':14}, title='tractography method')
    plt.setp(legend.get_title(),fontsize=14, fontweight='bold')
    plt.tight_layout()

    plt.show()



if __name__ == '__main__':
    dataset = 'ukb'     #dMRI data only for UKB
    which = 'total_length'    #avg_degree, total_length or total_density

    output_qball = 'ukb_totals.csv'
    output_dti = 'ukb_totals_dti.csv'
    df_qball = pd.read_csv(output_qball)
    df_dti = pd.read_csv(output_dti)


    phenotypes = pd.read_csv('../../ukb/phenotypes.csv')

    df_qball = df_qball.merge(phenotypes,on='eid')
    df_dti = df_dti.merge(phenotypes,on='eid')


    plotout(df_qball['age'], df_qball[which], df_dti[which],which)
