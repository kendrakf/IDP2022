#!/usr/bin/env python


''' Plot the colorful histograms for each amino acid in the folders
     '''

import sys
import glob
import argparse
from os.path import exists
from tqdm import tqdm
from Scripts.plot_1D_hist_w_regions import main as hist_parceled_cs

def main():
  
  folders = glob.glob("011_*/*_2k.out")
  for i in folders:
    aa = i.split("/")[-1][0]
    if aa != "G":
      #print(aa)
      hist_parceled_cs('{}'.format(i), 
                       ('015_1D_histogram_plots/'
                       'fig_{}.pdf'.format(aa)), 
                       '{}'.format(aa),
                       'frederick',
                       'redblue',
                       1.00,
                      )
 

if __name__ == "__main__":

  main()

