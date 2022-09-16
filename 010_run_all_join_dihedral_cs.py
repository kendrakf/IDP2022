#!/usr/bin/env python


''' Join the information from the folders with PDB files and files containing PPM predictions. 
    The result is a single file that contains the phi, psi, and all CS for the analyzed amino acid. 
     '''

import sys
import glob
import argparse
from Scripts.join_dihedral_cs import main as join_dihedral_cs
from Bio.Data.IUPACData import protein_letters_1to3 as one2three

def main():

  for i in list("ADEFGHIKLMNPQRSTVYWXZ"): # amino acid for each peptide
                                          # except CYS for which there is no prediction
    # prepare filenames and input options
    input_folder_pdb = "005_{}_fm2_fixbb".format(i)
    input_folder_cs = "008_{}_fm2_cs".format(i)
    output_path = "011_{}_rama_cs/{}_2k.out".format(i, i)
    if i == "X": # X is pre-proline: A followed by P 
      aa1 = "A"
      aa3 = "ALA"
      print(output_path, aa1, aa3)
    elif i == "Z": # Z is proline followed by proline
      aa1 = "P"
      aa3 = "PRO"
      print(output_path, aa1, aa3)
    else: # for all other regular amino acids
      aa1 = i
      aa3 = '{}'.format(one2three[aa1].upper())
      print(output_path, aa1, aa3)

    # run the main script for the central AA in 
    #   the A(25)-X-A(25) peptide ensemble
    join_dihedral_cs(
                  input_folder_pdb,
                  input_folder_cs,
                  output_path,
                  26,
                  26,
                  aa3,
                  )

if __name__ == "__main__":
  
  main()

