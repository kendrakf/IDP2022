#!/usr/bin/env python


''' runs the two/three scripts that take in the in chemical shift of the ensemble,
    add gaussian linewidth, and sum up everything to result in a simulated spectrum.
    '''

import os
import sys
import glob
import argparse
from Scripts.cs_in_spectrum_out import main as cs_in_spectrum_out
from Scripts.spectrum_in_line_out import main as spectrum_in_line_out


def main():

  output_folder = '017_gaussian_simulations'
  print(output_folder)

  if not os.path.exists(output_folder):
      os.makedirs(output_folder)

  ### First write out a two column file: 1) ppm 2) peak height ###
  ### Output file has the ending ".spectrum" ###
  cs_in_spectrum_out('011_T_rama_cs/T_2k.out', output_folder, 'T', 1.900)
  #cs_in_spectrum_out('011_T_rama_cs/T_2k.out', output_folder, 'T', 2.355)
  #cs_in_spectrum_out('011_T_rama_cs/T_2k.out', output_folder, 'T', 3.000)

  cs_in_spectrum_out('011_V_rama_cs/V_2k.out', output_folder, 'V', 1.900)
  #cs_in_spectrum_out('011_V_rama_cs/V_2k.out', output_folder, 'V', 2.355)
  #cs_in_spectrum_out('011_V_rama_cs/V_2k.out', output_folder, 'V', 3.000)

  ### Use the ".spectrum" files to draw the single-line plots ###
  spectrum_in_line_out('017_gaussian_simulations', '017_gaussian_simulations', 'long')
  spectrum_in_line_out('017_gaussian_simulations', '017_gaussian_simulations', 'long')


if __name__ == "__main__":

  main()

