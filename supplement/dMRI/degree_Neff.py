"""
Plot individually fitted Neff as a function of 
average degree.

"""


import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import spearmanr, binned_statistic



if __name__ == '__main__':
    dataset = 'ukb' #only available for ukb cuz only dataset with dMRI
    which = 'avg_degree'    #avg_degree, length or density

    output_qball = 'ukb_totals.csv'
    df_qball = pd.read_csv(output_qball)

    Neffs = pd.read_csv('ukb_Neffs.csv')


    data = df_qball.merge(Neffs,on='eid')
    rho, pval = spearmanr(data['avg_degree'], data['Neff'])
    plt.scatter(data['avg_degree'], data['Neff'], label='$\\rho=$ {0:.2f} ({1:.2E})\nN = {2}'.format(rho, pval, len(data['avg_degree'])), color='m')

        
    dataset_title = 'UK Biobank'
    plt.title(dataset_title, fontsize=18, fontweight='bold')
    plt.xlabel('average degree', fontsize=18)


    plt.ylabel('$N_{\mathrm{eff}}$', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(prop={'size':14})
    plt.tight_layout()

    plt.show()
