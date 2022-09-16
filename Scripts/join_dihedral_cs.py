#!/usr/bin/env python


''' Input is a folder of PDB files.
    Output is a list of all dihedrals '''

from __future__ import print_function
from __future__ import with_statement
from __future__ import division

__author__ = 'Jaka Kragelj'
__email__  = 'jaka.kragelj@utsouthwestern.edu'
__license__ = 'Creative Commons Attribute By - http://creativecommons.org/licenses/by/3.0/us/'''

import os
import sys
import csv
import glob
import math
from tqdm import tqdm
import pandas as pd
import argparse
import glob
from Bio.PDB import *


def main(pdb_folder, cs_folder, output_filename, begin, end, resi_choice):
  
  ### some safety precautions ###
  resi_choice = resi_choice.upper()
  ### get the *.pdb file names ###
  pdbs = glob.glob(pdb_folder + "/*.pdb")
  parser = PDBParser(QUIET=True)
  ### THE DATAFRAME ###
  columns = ['#res','aa','phi','psi','ca', 'cb', 'c', 'n']
  df = pd.DataFrame(columns=columns)
  ### process each PDB separately ###
  for pdb in pdbs: # UNCOMMENT when running with sbatch 
# for pdb in tqdm(pdbs): # COMMENT IN when running sbatch 
                         # (shows progress in the terminal)
    filename_core = pdb.split("/")[-1].split(".")[0]
    cs_path = cs_folder + "/" + filename_core + ".cs"
    ### check the *.cs file name ###
    if not os.path.isfile(cs_path):
      print("Error for file {}".format(filename_core))
      print("The PDB file and CS file don't have matching file cores\n")
      continue

    ### get the CS data ###
    current_aas = dict()
    CAs=dict(); CBs=dict(); COs=dict(); Ns=dict()
    chem_shifts = {'CA':CAs, 'CB':CBs, 'CO':COs, 'N':Ns}
    with open(cs_path, 'r') as handle:
      body_flag = False
      for line in handle:
        if body_flag == True:
          if line.strip() != "stop_": # stop_ means the CS data has ended
            atm_nr, res_nr, aa, at_id, at_id2, cs, dot, one = line.strip().split()
            if at_id == 'CA':
              res_nr = int(res_nr)
              chem_shifts['CA'][res_nr]=float(cs)
              current_aas[res_nr]=aa
            elif at_id == 'CB':
              res_nr = int(res_nr)
              chem_shifts['CB'][res_nr]=float(cs)
              current_aas[res_nr]=aa
            elif at_id == 'C': # "C" is used for CO/C' atom!
              res_nr = int(res_nr)
              chem_shifts['CO'][res_nr]=float(cs)
              current_aas[res_nr]=aa
            elif at_id == 'N':
              res_nr = int(res_nr)
              chem_shifts['N'][res_nr]=float(cs)
              current_aas[res_nr]=aa
        elif line.strip() == "_Chem_shift_ambiguity_code": # start signal
          body_flag = True
#   print(chem_shifts)

    ### some corrections to input ###
    ### This code block is perhaps not necessary. ###
    if begin == 1:
      print("Can't do join for the first amino"
            " acid because usually chemical shifts"
            " are not predicted for the first AA."
            " Moving starting point to res 2.")
      begin = 2
    if end == max(current_aas.keys()):
      print("Can't do join for the last amino"
            " acid because usually chemical shifts"
            " are not predicted for the first AA."
            " Moving starting point to res max().")
      end = max(current_aas.keys())
    if begin < min(current_aas.keys()):
      begin = min(current_aas.keys())
    if end > max(current_aas.keys()):
        end = max(current_aas.keys())
    
    ### get the dihedral data ...      ###
    ### ... combine with CS data ...   ###
    ### ... and store in a dictionary. ###
    data = parser.get_structure('conformer', pdb)
    for model in data.get_models():
      for chain in model.get_chains():
        res_names = []
        for residue in chain.get_residues():
            res_names.append(residue.get_resname())
        poly = Polypeptide.Polypeptide(chain)
        phis_psis = poly.get_phi_psi_list()
        res_list = list(chain.get_residues())
        for m, n in zip(res_list,phis_psis):
          i = int(m.id[1])
          if begin <= i <= end: # syniclein: 1-140
            if m.get_resname() == resi_choice:
#             print(i, m.get_resname(),n[0]*180/math.pi, n[1]*180/math.pi, current_aas[i], chem_shifts['CA'][i])
              if m.get_resname() == "GLY":
                df.loc[len(df)] = [
                                   i, 
                                   m.get_resname(),
                                   n[0]*180/math.pi, 
                                   n[1]*180/math.pi, 
                                   chem_shifts['CA'][i],
                                   0.000,
                                   chem_shifts['CO'][i],
                                   chem_shifts['N'][i],
                                  ]
              elif m.get_resname() == "PRO":
                df.loc[len(df)] = [
                                   i, 
                                   m.get_resname(),
                                   n[0]*180/math.pi, 
                                   n[1]*180/math.pi, 
                                   chem_shifts['CA'][i],
                                   chem_shifts['CB'][i],
                                   chem_shifts['CO'][i],
                                   0.000,
                                  ]
              else:
                df.loc[len(df)] = [
                                   i, 
                                   m.get_resname(),
                                   n[0]*180/math.pi, 
                                   n[1]*180/math.pi, 
                                   chem_shifts['CA'][i],
                                   chem_shifts['CB'][i],
                                   chem_shifts['CO'][i],
                                   chem_shifts['N'][i],
                                  ]


  ### write to file ###
  df.to_csv(output_filename, 
            sep="\t",
            index=False,
            float_format='%8.3f',
            quoting=csv.QUOTE_NONE,
            escapechar="\\"
            )


if __name__ == "__main__":

  parser = argparse.ArgumentParser(description="take in .*pdb and *.cs and write to file.")
  parser.add_argument("input_pdb_folder", help="input folder containing PDBs")
  parser.add_argument("input_cs_folder", help="input folder containing files of same name ending with .cs")
  parser.add_argument("output_filename", help="output file to write to.")
# parser.add_argument("-a","--atom", default="ALL", help="Specify: ALL. Not implemented: CA, CB, CO, N")
  parser.add_argument("-b","--begin_number", type=int, default=2, help="first residue number")
  parser.add_argument("-e","--end_number", type=int, default=138, help="first residue number")
  parser.add_argument("-r","--residue_type", type=str, default='THR', help="Uppercase amino acid")
  args = parser.parse_args()

  main(args.input_pdb_folder, 
       args.input_cs_folder,
       args.output_filename,
       args.begin_number,
       args.end_number,
       args.residue_type,
       )
