#!/usr/bin/env python


''' Plot histograms of the chemical shifts for a single atom and
    color them according to the regions of the ramachandran plot. '''

import sys
import glob
import math
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as pltticker
from matplotlib import colors
import seaborn as sns

 
#################################################
### GENERAL PROPERTIES x-limits, rama regions ###
#################################################

na=np.nan
atoms = ['ca','cb','c']
aa_spec_ranges={
                # Specify the x-axis limits.
                # Ca & Cb adjusted to 17ppm to match Fig5
                #aa1 Ca      Cb      C
                'A':[[44,61],[11,28],[170,184]],
                'R':[[47,64],[22,39],[169,183]],
                'N':[[44,61],[30,47],[168,182]],
                'D':[[45,62],[33,50],[169,183]],
#               'C':[[50,67],[22,39],[168,182]],
                'Q':[[47,64],[22,39],[169,183]],
                'E':[[49,66],[22,39],[169,183]],
                'G':[[36,53],[na,na],[167,181]],
                'H':[[49,66],[23,40],[168,182]],
                'I':[[54,71],[31,48],[168,182]],
                'L':[[45,62],[33,50],[170,184]],
                'K':[[49,66],[26,43],[169,183]],
                'M':[[47,64],[26,43],[169,183]],
                'F':[[50,67],[31,48],[168,182]],
                'P':[[57,74],[25,42],[170,184]],
                'S':[[50,67],[55,72],[167,181]],
                'T':[[54,71],[60,77],[168,182]],
                'W':[[48,65],[21,38],[169,183]],
                'Y':[[50,67],[30,47],[168,182]],
                'V':[[54,71],[25,42],[169,183]],
                }

### specification of how the ramachandran   ###
### space is separated into parcels/regions ###
rama_definitions = {
  ### selection: Frederick lab ###
  # ralpha and rturnI can overlap since
  # there is statement later on to exclude 
  # ralpha from rtunrI. This exclusion doesn't 
  # affect the "blackledge: and "bax" selection
  # since there is no overlap in the later two. 
  "frederick":
  {
  'rbetaA':(-180.0, -90.0, 105.0, 180.0),
  'rbetaB':(-180.0,-105.0,-180.0,-135.0),
  'rppiiA':( -90.0,   0.0, 105.0, 180.0),
  'rppiiB':(-105.0,   0.0,-180.0,-135.0),
  'rdelta':(-180.0,   0.0,  45.0, 105.0),
  'rlhelx':(   0.0, 180.0,-180.0, 180.0),
  'rturnI':(-180.0,   0.0,-135.0,  45.0), # turn is wider alpha
  'ralpha':( -82.5, -37.5, -60.0, -20.0), # alph is strictly alpha
  },
  ### selection: Blackledge lab ###
  "blackledge":
  {
  'rbetaA':(-180.0, -90.0,  45.0, 180.0),
  'rbetaB':(-180.0, -90.0,-180.0,-120.0),
  'rppiiA':( -90.0,   0.0,  45.0, 180.0),
  'rppiiB':( -90.0,   0.0,-180.0,-120.0),
  'rdelta':(   0.0,   0.0,   0.0,   0.0),
  'rlhelx':(   0.0, 180.0,-180.0, 180.0),
  'rturnI':(   0.0,   0.0,   0.0,   0.0),
  'ralpha':(-180.0,   0.0,-120.0,  45.0),
  },
  ### selection: Bax lab ###
  "bax":
  {
  'rbetaA':(-180, -90,  90, 180),
  'rbetaB':(   0,   0,   0,   0),
  'rppiiA':( -90, -45, 105, 180),
  'rppiiB':(   0,   0,   0,   0),
  'rdelta':(   0,   0,   0,   0),
  'rlhelx':(  45,  75,  15,  60),
  'rturnI':(-135, -75, -15,  30),
  'ralpha':( -90, -45, -60, -15),
  }}

# the "rama_definitions" and "label_definitions" can be extended
# to include more regions (e.g. separate epsilono and delta'). 
labels_definitions = {
        'frederick': ['beta','ppii','ddgg','left','turn','alph'],
        'blackledge':['beta','ppii','ddgg','left','turn','alph'],
        'bax':       ['beta','ppii','ddgg','left','turn','alph'],
        }

# Order of the colors: 'beta','ppii','ddgg','left','turn','alph'
seaborn_colorblind = [colors.to_hex(i)for i in sns.color_palette('colorblind',6)]
seaborn_husl = [colors.to_hex(i)for i in sns.color_palette('husl',6)]
color_definitions = {
       'clown':  ['#0000f2','#00a88d','#ebb200','#e88c0b','#f0601e','#f00000'],
       'viridis':['#440356','#3e4989','#297b8e','#23a983','#7fd34e','#fbe723'],
       'plasma': ['#0d0887','#6900a8','#b02a90','#e06363','#fca537','#f0f724'],
       'redblue':['#1b1b64', '#5e93dc', '#b3dfff', '#f3ecb7', '#eeaa38', '#a60000'],
       'seaborn': seaborn_colorblind,
       'husl': seaborn_husl,
       }


################################################
###     main() function: data processing     ###
################################################

def main(input_file_path, output_file_path, aa, definition, color_palette, alpha_choice):


  ### read in regions ###
  df = pd.read_csv(input_file_path, sep="\t", header=0)
  # the columns/header of the input should be:
  # res  aa  phi  psi ca  cb  c  n

  ### order of ramachandran elements ###
  labels=labels_definitions[definition]

  ### which definition of regions in Ramachandran  ###
  ### space have you selected while providing args ###
  rchoice = rama_definitions[definition]

  ### choice of the color palette ###
  col = color_definitions[color_palette]
  
  ### extract plot subregions ###
  # "r" stands for ramachandran
  # regions: beta-strand, polyproline helix ii
  #          beta-turn I, alpha-helix, 
  #          left handed helix (all phi > 0)
  # !!! top and left corner are always >= or <=
  #     unless they reach the bottom/right axis
  #     (i.e. betaB, ppiiB, left)
  data_beta = pd.concat([
              df.loc[
                     (df['phi'] >= rchoice['rbetaA'][0]) & 
                     (df['phi'] <  rchoice['rbetaA'][1]) & 
                     (df['psi'] >  rchoice['rbetaA'][2]) & 
                     (df['psi'] <= rchoice['rbetaA'][3])
                     ],
              df.loc[
                     (df['phi'] >= rchoice['rbetaB'][0]) & 
                     (df['phi'] <  rchoice['rbetaB'][1]) & 
                     (df['psi'] >= rchoice['rbetaB'][2]) & 
                     (df['psi'] <= rchoice['rbetaB'][3])
                     ]])
  data_ppii = pd.concat([
              df.loc[
                     (df['phi'] >= rchoice['rppiiA'][0]) & 
                     (df['phi'] <  rchoice['rppiiA'][1]) & 
                     (df['psi'] >  rchoice['rppiiA'][2]) & 
                     (df['psi'] <= rchoice['rppiiA'][3])
                     ],
              df.loc[
                     (df['phi'] >= rchoice['rppiiB'][0]) & 
                     (df['phi'] <  rchoice['rppiiB'][1]) & 
                     (df['psi'] >= rchoice['rppiiB'][2]) & 
                     (df['psi'] <= rchoice['rppiiB'][3])
                     ]])
  data_ddgg = df.loc[
                     (df['phi'] >= rchoice['rdelta'][0]) & 
                     (df['phi'] <  rchoice['rdelta'][1]) & 
                     (df['psi'] >  rchoice['rdelta'][2]) & 
                     (df['psi'] <= rchoice['rdelta'][3])
                     ]
  data_turn = df.loc[
                    (
                     (df['phi'] >= rchoice['rturnI'][0]) & 
                     (df['phi'] <  rchoice['ralpha'][0]) & 
                     (df['psi'] >= rchoice['rturnI'][2]) & 
                     (df['psi'] <= rchoice['rturnI'][3])  
                    ) | (
                     (df['phi'] >= rchoice['ralpha'][1]) & 
                     (df['phi'] <  rchoice['rturnI'][1]) & 
                     (df['psi'] >= rchoice['rturnI'][2]) & 
                     (df['psi'] <= rchoice['rturnI'][3])  
                    ) | (
                     (df['phi'] >= rchoice['ralpha'][0]) & 
                     (df['phi'] <  rchoice['ralpha'][1]) & 
                     (df['psi'] >= rchoice['rturnI'][2]) & 
                     (df['psi'] <= rchoice['ralpha'][2])   
                    ) | (
                     (df['phi'] >= rchoice['ralpha'][0]) & 
                     (df['phi'] <  rchoice['ralpha'][1]) & 
                     (df['psi'] >  rchoice['ralpha'][3]) &
                     (df['psi'] <= rchoice['rturnI'][3])  
                    )
                    ]
  data_alph = df.loc[
                     (df['phi'] >= rchoice['ralpha'][0]) & 
                     (df['phi'] <  rchoice['ralpha'][1]) & 
                     (df['psi'] >  rchoice['ralpha'][2]) & 
                     (df['psi'] <= rchoice['ralpha'][3])
                     ]
  data_left = df.loc[
                     (df['phi'] >= rchoice['rlhelx'][0]) & 
                     (df['phi'] <= rchoice['rlhelx'][1]) & 
                     (df['psi'] >= rchoice['rlhelx'][2]) & 
                     (df['psi'] <= rchoice['rlhelx'][3])
                     ]

  # since DataFrame doesn't work for different lengths 
  # of series for CS of each region: 
  #    reorganize into atom specific dict with columns
  #    representing cs predictions for the regions
  cs_ca = ({'beta': data_beta['ca'].values.tolist(),
            'ppii': data_ppii['ca'].values.tolist(),
            'ddgg': data_ddgg['ca'].values.tolist(),
            'left': data_left['ca'].values.tolist(),
            'turn': data_turn['ca'].values.tolist(),
            'alph': data_alph['ca'].values.tolist(),
            })
  cs_cb = ({'beta': data_beta['cb'].values.tolist(),
            'ppii': data_ppii['cb'].values.tolist(),
            'ddgg': data_ddgg['cb'].values.tolist(),
            'left': data_left['cb'].values.tolist(),
            'turn': data_turn['cb'].values.tolist(),
            'alph': data_alph['cb'].values.tolist(),
            })
  cs_co = ({'beta': data_beta['c'].values.tolist(),
            'ppii': data_ppii['c'].values.tolist(),
            'ddgg': data_ddgg['c'].values.tolist(),
            'left': data_left['c'].values.tolist(),
            'turn': data_turn['c'].values.tolist(),
            'alph': data_alph['c'].values.tolist(),
            })

  cs_all = {'ca':cs_ca, 
            'cb':cs_cb, 
            'c' :cs_co,
            }

# print("\nbeta"); print(data_beta); print("\nppii"); print(data_ppii); print("\ndelta-gamma"); print(data_ddgg); print("\nturn"); print(data_turn); print("\nalph"); print(data_alph); print("\nleft"); print(data_left); sys.exit()

  # CALCULATION OF HISTOGRAM AND PLOTTING #
  plot_types = [False,True] # Histograms are either True/Stacked or False/Overlapped
  for plot_type in plot_types:
    
    # initiate plotting
    fig, axs = plt.subplots(1, 3, sharey=False, tight_layout=True,figsize=(9.75,1.50))

    # specify bind widths
    min_max=aa_spec_ranges[aa]
    bin_width=0.25
    bins_array=[np.arange(i[0],i[1]+1,bin_width) for i in min_max]

    # PLOTTING all three atoms: CA, CB, C
    for i in range(len(atoms)): 
      # HIST PLOT
      # hist() uses plot_type variable to determine if stacked or not
      h, bins, patches = axs[i].hist(list(cs_all[atoms[i]].values()),
                            bins=bins_array[i],
                            density=False, 
                            histtype='stepfilled', 
                            stacked=plot_type,
                            color=col,
                            alpha=alpha_choice,
                            label=labels)

      stacked_data = [j for sub in list(cs_all[atoms[i]].values()) for j in sub]
      axs[i].hist(stacked_data, # SUM PLOT: line histogram plot
                  bins=bins_array[i], 
                  density=False, 
                  histtype='step', 
                  stacked=False,
                  color='black',
                  alpha=1.0)

      ### FORMATTING ###
      axs[i].set_xlim(min_max[i])
      mtick_spacing=5
      axs[i].xaxis.set_major_locator(pltticker.MultipleLocator(mtick_spacing))
      axs[i].xaxis.set_minor_locator(pltticker.MultipleLocator(1))
      axs[i].yaxis.set_ticks([])

      # other formatting
      axs[i].invert_xaxis()
      for line in ['top','left','right']:
        axs[i].spines[line].set_visible(False)
      for line in ['bottom']:
        axs[i].spines[line].set_visible(True)

    # SAVE PLOT: save figure with _stacked suffix if stacked
    if not plot_type:
      print("\nsaved: {}".format(output_file_path))
      plt.savefig(output_file_path)
    else:
      print("saved: {}\n".format(output_file_path.split('.')[0] + "_stacked.pdf"))
      plt.savefig(output_file_path.split('.')[0] + "_stacked.pdf")
    
    plt.clf()
  plt.close()


if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='this is a template script that prints our file contents.')
  parser.add_argument('input_filename', 
                      type=str,
                      help='input file'
                      )
  parser.add_argument('output_filename', 
                      type=str,
                      help='output file to write to.'
                      )
  parser.add_argument('-a',
                      '--amino_acid', 
                      default='T',
                      type=str,
                      help='Which amino acid? T or V?'
                      )
  parser.add_argument('-r',
                      '--ramachandran_definition', 
                      default='frederick',
                      type=str,
                      help=('Three options: frederick, blackledge, '
                            'bax. Without specifying "frederick" is '
                            'used as a default.')
                      )
  parser.add_argument('-c',
                      '--color_palette', 
                      default='viridis',
                      type=str,
                      help=('Three options: viridis, plasma, clown.\n'
                            '"viridis" used as a default.')
                      )
  parser.add_argument('-t',
                      '--transparency', 
                      default=0.8,
                      type=float,
                      help=('Transparency/alpha: 0.0 - 1.0')
                      )
  args = parser.parse_args()

  main(args.input_filename, args.output_filename, args.amino_acid, args.ramachandran_definition, args.color_palette, args.transparency)
