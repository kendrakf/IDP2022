#!/usr/bin/env python

''' Description: takes in spectrum coordinates from a folder
                 returns line graph output in png form and pdf form'''


import sys
import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import glob

aa_spec_ranges={
                     #CA/CB  #CA/CB  #CO       #N
                'T':[[54,71],[60,77],[168,182],[102,123]],
                'V':[[54,71],[25,42],[169,183],[108,129]],
                }
atom_reference = {'CAs':0,'CBs':1,'COs':2,'Ns':3}

def main(input_folder, output_folder, shape):

    spectra_files = glob.glob("{}/*.spectrum".format(input_folder))

    for input_path in spectra_files:
      print("Plotting {}".format(input_path))
      # prepare
      filename = input_path.split("/")[-1]
      aa, atoms, tmp1 = filename.split("_")
      min_max = aa_spec_ranges[aa][atom_reference[atoms]]
      # go
      ax = plt.axes()
      if shape == 'square':
          plt.gcf().set_size_inches(3.5, 3.5)
      if shape == 'long':
          plt.gcf().set_size_inches(7, 3.5)
      ax.spines['bottom'].set_visible(True)
      ax.spines['top'].set_visible(False)
      ax.spines['right'].set_visible(False)
      ax.spines['left'].set_visible(False)
      ax.autoscale(axis='x', tight=True)
      mpl.rc('font', family='Arial')

      df = pd.read_csv(input_path, delim_whitespace=True, header=None, names=['xvalues', 'yvalues'])

      ax.plot('xvalues', 'yvalues', data=df, linewidth=1)

      xlocs, xlabels = plt.xticks()
      plt.xticks(ticks=xlocs, fontfamily='Arial', fontsize=10)
      plt.yticks([])
      
      # new section
      tick_spacing=5
      tick_min = (min_max[0]+5)//5*5
      tick_max = (min_max[1]+0)//5*5
      min_max_ticks=np.arange(tick_min, tick_max + 1,tick_spacing)
      ax.set_xlim(min_max)
      ax.xaxis.set_ticks(min_max_ticks)
      ax.xaxis.set_minor_locator(plticker.MultipleLocator(1))
      ax.invert_xaxis()

      input_filename = input_path.split("/")[-1]
      output_name_core = input_filename.split(".")[0]
      plt.savefig(output_folder + "/figure_" + output_name_core + ".pdf")
      plt.savefig(output_folder + "/figure_" + output_name_core + ".png", dpi=360)
      plt.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='outputs coordinates in txt file as line graph (both pdf and png).')
    parser.add_argument(
                        '-i',
                        '--input_folder',
                        type=str,
                        required=True,
                        help="input text file containing columns: x-value and y-value.")
    parser.add_argument(
                        '-o',
                        '--output_folder',
                        type=str,
                        required=True,
                        help='output path to write to (do not specify type).')
    parser.add_argument(
                        '-s',
                        '--shape',
                        type=str,
                        required=False,
                        default='long',
                        help='output graph: \'square\' for 3.5 x 3.5 in or \'long\' for 7 x 3.5 in (default=long).')
    args = parser.parse_args()

    main(args.input_folder, args.output_folder, args.shape)
