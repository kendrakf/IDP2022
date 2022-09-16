#!/usr/bin/env python


''' Take in file with eight columns (#nr #THR #phi #psi #CA #CB #CO #N) 
  and simulate a 1D for all atoms in the input.
  You can select a region of the protein. '''

import os
import sys
import glob
import argparse
import numpy as np
from tqdm import tqdm


def gaussian(x, mu, sig):
  return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)

aa_spec_ranges={
                'T':{'CAs':[51,73],'CBs':[58,80],'COs':[163,185],'Ns':[104,126]},
                'V':{'CAs':[55,69],'CBs':[25,39],'COs':[167,181],'Ns':[110,132]},
                }


def main(input_file_path, output_folder, aa, peak_width):

  if os.path.isdir(input_file_path):
    print("provide a path to a file, not a folder.")
    sys.exit()

  CAs=[]; CBs=[]; COs=[]; Ns=[]

  with open(input_file_path, "r") as handle:
    for line in handle:
      if "#" not in line:
        res_nr, res_id, phi, psi, tmp_ca, tmp_cb, tmp_co, tmp_n = line.strip().split()
        CAs.append(float(tmp_ca))
        CBs.append(float(tmp_cb))
        COs.append(float(tmp_co))
        Ns.append(float(tmp_n))

  # Store input in a dictionary
  all_shifts = {"CAs":CAs, "CBs":CBs, "COs":COs, "Ns":Ns}
  # Store partial spectra after applying gaussian
  all_simulated_spectra = {"CAs":[], "CBs":[], "COs":[], "Ns":[]}
  # Store sum of the partial spectra - the result
  sum_spectra = {}


  for atom_type in all_shifts.keys():
    print("Preparing {}".format(atom_type))

    # x-values
    x1, xn = aa_spec_ranges[aa][atom_type]
    x_values = np.arange(xn, x1, -0.01) # minus 0.01!!!
  
    # sigma
    # [!] peak_width is now specified in the function definition
    # peak_width = 2.355 #ppm
    peak_sigma = peak_width / 2.355
    sigma_cb = np.array([peak_sigma]*len(x_values))

    for cs in tqdm(all_shifts[atom_type]):
      cs_array = np.array([cs]*len(x_values))
      sigma_array = np.array([peak_sigma]*len(x_values))
      partial_spectrum = gaussian(x_values, cs_array, sigma_array)
      all_simulated_spectra[atom_type].append(partial_spectrum)
    
# print("\nSumming and writting")
# for atom_type in all_shifts.keys():
    tmp_np_sim_spectra = np.array(all_simulated_spectra[atom_type])
    sum_spectra[atom_type] = tmp_np_sim_spectra.sum(axis=0)


# for atom_type in all_shifts.keys():
    atom_output_file_path = "{}/{}_{}_{}ppm.spectrum".format(
                             output_folder,aa,atom_type, str(peak_width).replace(".","p"))
    np.savetxt(atom_output_file_path, np.array([x_values,sum_spectra[atom_type]]).transpose())
  


if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Read a folder with SPARTA+ pred files, get CB shifts, and simulate gaussians.')
  parser.add_argument('-i',
            '--input_file_path',
            type=str,
            required=True,
            help='Folder containing many .cs files')
  parser.add_argument(
            '-o',
            '--output_folder',
            type=str,
            required=True,
            help='Output folder to write to.'
            )
  parser.add_argument(
            '-r',
            '--aa',
            type=str,
            default='T',
            help='single letter aa code'
            )
  parser.add_argument(
            '-w',
            '--full_width_half_max',
            type=float,
            default=2.355,
            help='Default FWHM is 2.355 ppm (sigma = 1ppm).')
  args = parser.parse_args()

  main(args.input_file_path, args.output_folder, args.aa, args.full_width_half_max)
