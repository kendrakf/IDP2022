#!/usr/bin/env python


''' plot DARR-like/TEDOR-like 2D 
    histogram of the chemical shifts '''

import sys
import glob
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as pltticker
from matplotlib import colors
from statistics import mean
from Bio.Data.IUPACData import protein_letters_1to3 as one2three

atom_key={
          "ca":0,
          "cb":1,
          "co":2,
          "n":3,
         }

nonCys_nonGly_aas={"A", "D", "E", "F", "H", "I", "K", "L", "M", 
                   "N", "P", "Q", "R", "S", "T", "V", "Y", "W"}

nonCys_aas={"A", "D", "E", "F", "G", "H", "I", "K", "L", "M", 
            "N", "P", "Q", "R", "S", "T", "V", "Y", "W"}

na=np.nan

aa_spec_ranges={
                #aa1 Ca      Cb      C         N
                # Ca & Cb adjusted to 17ppm to match Fig5
                'A':[[44,61],[11,28],[170,184],[112,133]],
                'R':[[47,64],[22,39],[169,183],[112,133]],
                'N':[[44,61],[30,47],[168,182],[106,127]],
                'D':[[45,62],[33,50],[169,183],[109,130]],
#               'C':[[50,67],[22,39],[168,182],[105,126]],
                'Q':[[47,64],[22,39],[169,183],[108,129]],
                'E':[[49,66],[22,39],[169,183],[109,130]],
                'G':[[36,53],[na,na],[167,181],[ 98,119]],
                'H':[[49,66],[23,40],[168,182],[107,128]],
                'I':[[54,71],[31,48],[168,182],[109,130]],
                'L':[[45,62],[33,50],[170,184],[110,131]],
                'K':[[49,66],[26,43],[169,183],[109,130]],
                'M':[[47,64],[26,43],[169,183],[108,129]],
                'F':[[50,67],[31,48],[168,182],[107,128]],
                'P':[[57,74],[25,42],[170,184],[ na, na]],
                'S':[[50,67],[55,72],[167,181],[104,125]],
                'T':[[54,71],[60,77],[168,182],[102,123]],
                'W':[[48,65],[21,38],[169,183],[109,130]],
                'Y':[[50,67],[30,47],[168,182],[108,129]],
                'V':[[54,71],[25,42],[169,183],[108,129]],
                }


aa_atom_groups={
                'A':[['co','ca'],['ca','cb']],
#               'C':[['co','ca'],['ca','cb']],
                'D':[['co','ca'],['ca','cb']],
                'E':[['co','ca'],['ca','cb']],
                'F':[['co','ca'],['ca','cb']],
                'G':[['co','ca']],
                'H':[['co','ca'],['ca','cb']],
                'I':[['co','ca'],['ca','cb']],
                'K':[['co','ca'],['ca','cb']],
                'L':[['co','ca'],['ca','cb']],
                'M':[['co','ca'],['ca','cb']],
                'N':[['co','ca'],['ca','cb']],
                'P':[['co','ca'],['ca','cb']],
                'Q':[['co','ca'],['ca','cb']],
                'R':[['co','ca'],['ca','cb']],
                'S':[['co','ca'],['ca','cb']],
                'T':[['co','ca'],['ca','cb']],
                'V':[['co','ca'],['ca','cb']],
                'Y':[['co','ca'],['ca','cb']],
                'W':[['co','ca'],['ca','cb']],
               }

def main(atom_a, atom_b, above_diag, output_filename):
  #####################################
  ### SECTION BELOW CAN BE MODIFIED ###
  # find input files, name output - can be modified #
  input_file_paths=sorted(glob.glob("011_*_rama_cs/*2k.out"))
  # prepare for plotting layout and figure size #
  fig_grid_shape=(6,4) # modify for different grid arrangement of subplots
  fig, axs = plt.subplots(fig_grid_shape[0],fig_grid_shape[1],
                          figsize=(7,9.5), tight_layout=True)
                          # figsize is in inches, modify for different size

  # begin loop that goes over each aa/input/subplot #
  for subplot_counter, input_file in enumerate(input_file_paths):
    # also, extract amino acid information #
    # get the first letter of the input file - modify manually if necessary #
    aaX = input_file.split("/")[-1][0] # modify this to get aa1 info

    ###################################
    ### AUTOMATIC PART FROM HERE ON ###
    # skip G and P when the specified atom/CS is not there
    atom_a = atom_a.lower(); atom_b = atom_b.lower()
    if aaX == "G" and (atom_a == 'cb' or atom_b == 'cb'):
      print(subplot_counter, aaX, input_file, "skipping Gly: no Cbeta")
    elif aaX == "P" and (atom_a == 'n' or atom_b == 'n'):
      print(subplot_counter, aaX, input_file, "skipping Pro: no N CS information")
    elif aaX == "Z" and (atom_a == 'n' or atom_b == 'n'):
      print(subplot_counter, aaX, input_file, "skipping Pro-Pro: no N CS information")
    else:
      # assign the right aa1 to ala-pro (X) and pro-pro (Z)
      if aaX in nonCys_aas:
        aa1 = aaX
      elif aaX == "X":
        aa1 = "A"
      elif aaX == "Z":
        aa1 = "P"
      else:
        print("Error [{}]: unknown single letter amino acid code "
              "(or it was a Cys).".format(aaX))
      # read in CS for each file #
      print(subplot_counter, aaX, input_file)
      CAs=[]; CBs=[]; COs=[]; Ns=[] 
      with open(input_file, "r") as handle:
        for line in handle:
          if "#" not in line:
            if aa1 == 'G':
              num, aa3, phi_i, psi, ca_i, cb_i, co_i, n_i = line.strip().split()
              CAs.append(float(ca_i))
              COs.append(float(co_i))
              Ns.append(float(n_i))
            else:
              num, aa3, phi_i, psi, ca_i, cb_i, co_i, n_i = line.strip().split()
              CAs.append(float(ca_i))
              CBs.append(float(cb_i))
              COs.append(float(co_i))
              Ns.append(float(n_i))

      # join CS of each atom into one list #
      if aa1 in nonCys_nonGly_aas:
        all_CS = [CAs, CBs, COs, Ns]
      elif aa1 == 'G':
        all_CS = [CAs, CAs, COs, Ns]
      elif aa1 == 'C':
        print('Skipping cysteine since PPM has no prediction for Cys.')
      else:
        print('Error: amino acid \"{}\" not implemented'.format(aa))
        sys.exit()

      # set atom selection so that lower average(ppm) is first (x-axis) ...
      # ... which draws above peaks above the diagonal. Higher avgerage(ppm) ...
      # ... draws the peaks as peaks below the diagonal. 
      if atom_a == 'n' or atom_b == 'n': # for DARR (nitrogen is for TEDOR)
        if atom_a == 'n': # N is first axis
          x_key=atom_key[atom_a]
          y_key=atom_key[atom_b]
        else:
          x_key=atom_key[atom_b]
          y_key=atom_key[atom_a]
      if (atom_a == 'ca' or atom_a == 'cb') and (atom_b == 'ca' or atom_b == 'cb'):
        if aa1 != 'S' and aa1 != 'T': # T & S behave just the opposite
          if mean(all_CS[atom_key[atom_a]]) > mean(all_CS[atom_key[atom_b]]) and above_diag:
            x_key=atom_key[atom_a]
            y_key=atom_key[atom_b]
          else:
            x_key=atom_key[atom_b]
            y_key=atom_key[atom_a]
        else: # for T & S
          if mean(all_CS[atom_key[atom_a]]) > mean(all_CS[atom_key[atom_b]]):
            if above_diag:
              x_key=atom_key[atom_b]
              y_key=atom_key[atom_a]
            else:
              x_key=atom_key[atom_a]
              y_key=atom_key[atom_b]
          else:
            if above_diag:
              x_key=atom_key[atom_a]
              y_key=atom_key[atom_b]
            else:
              x_key=atom_key[atom_b]
              y_key=atom_key[atom_a]
      else:
        if mean(all_CS[atom_key[atom_a]]) > mean(all_CS[atom_key[atom_b]]) and above_diag:
          x_key=atom_key[atom_a]
          y_key=atom_key[atom_b]
        else:
          x_key=atom_key[atom_b]
          y_key=atom_key[atom_a]

      # CALCULATION OF HISTOGRAM AND PLOTTING #
      # set bin parameters #
      bin_width=0.50
      mtick_spacing=5
      x_lims=np.array(aa_spec_ranges[aa1][x_key])
      y_lims=np.array(aa_spec_ranges[aa1][y_key])
      x_edges=np.arange(x_lims[0],x_lims[1]+bin_width,bin_width)
      y_edges=np.arange(y_lims[0],y_lims[1]+bin_width,bin_width)
      hw=bin_width/2
      center_bins_x=np.arange(x_lims[0]+hw,x_lims[1]+hw,bin_width)
      center_bins_y=np.arange(y_lims[0]+hw,y_lims[1]+hw,bin_width)

      # calculate counts for bins #
      H, npxedges, npyedges = np.histogram2d(all_CS[y_key], all_CS[x_key],
                                         bins=(y_edges,x_edges))
      # set contouring parameters #
      max_cont = np.max(H)
      ensemble_size = len(all_CS[y_key])
      min_cont = ensemble_size/2000*(5*bin_width) # ensemble_size/2000*(5units*bin_width)
                                      # formula would differ for a structured ensemble
      step_cont = max_cont/12 # 12 contour lines in a plot
      levels=np.arange(min_cont,max_cont,step_cont)

      # draw the subplot #
      m = subplot_counter//fig_grid_shape[1]
      n = subplot_counter%fig_grid_shape[1]
      c = axs[m][n].contour(x_edges[:-1], y_edges[:-1], H, levels, cmap='Oranges_r',linewidths=0.5)
      
      # adjust and format axes #
      axs[m][n].set_aspect('equal')
      axs[m][n].set_xlim(x_lims[0],x_lims[1])
      axs[m][n].set_ylim(y_lims[0],y_lims[1])
      axs[m][n].xaxis.set_major_locator(pltticker.MultipleLocator(mtick_spacing))
      axs[m][n].yaxis.set_major_locator(pltticker.MultipleLocator(mtick_spacing))
      axs[m][n].xaxis.set_minor_locator(pltticker.MultipleLocator(1))
      axs[m][n].yaxis.set_minor_locator(pltticker.MultipleLocator(1))
      axs[m][n].set_xlabel('ppm', labelpad=1.0)
      axs[m][n].set_ylabel('ppm', labelpad=1.0)
      axs[m][n].invert_xaxis()
      axs[m][n].invert_yaxis()

      # assign titles showed inside the subplots
      if aaX in nonCys_aas:
        inside_title=one2three[aa1].upper()
      elif aaX == "X":
        inside_title="ALA-PRO"
      elif aaX == "Z":
        inside_title="PRO-PRO"
      else:
        print("Error [{}]: unknown single letter amino acid code "
              "(or it was a Cys).".format(aaX))
      axs[m][n].text(0.95, 0.85, "{}".format(inside_title), ha="right", transform=axs[m][n].transAxes)

  # save plot
  plt.savefig(output_filename)
  print("Plotted: {}".format(output_filename))




if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Simulate correlation between a and b atoms.'
                                               ' Default for C atoms is DARR arrangement.'
                                               ' TEDOR is always plotted with N on x-axis.')
  parser.add_argument(
                      '-a', 
                      '--atom_a', 
                      type=str,
                      required=True,
                      help='x-axes centered on selected atom. Choices: ca, cb, co, n'
                      )
  parser.add_argument(
                      '-b', 
                      '--atom_b', 
                      type=str,
                      required=True,
                      help='y-axes centered on selected atom. Choices: ca, cb, co, n'
                      )
  parser.add_argument(
                      '-c', 
                      '--above_diagonal', 
                      type=bool,
                      default=True,
                      help='Add -c option to plot the peaks below the diagonal.'
                      )
  parser.add_argument('-o',
                      '--output_filename',
                      type=str,
                      default='test.pdf',
                      help="Output file to write to. Default: test.pdf")
  args = parser.parse_args()

  main(args.atom_a, args.atom_b, args.above_diagonal, args.output_filename)
