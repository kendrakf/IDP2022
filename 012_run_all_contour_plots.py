#!/usr/bin/env python


''' Run the "plot_contour_histogram_multi.py" script
    with different preset arguments so that to cover
    all atom combinations in contour plots. 
     '''

import sys
import glob
import argparse
from Scripts.plot_contour_histogram_multi import main as plot_contour

def main():

  plot_contour('ca', 'cb', True,  '013_contour_plots/figure_DARR_CaCb_top.pdf')
# plot_contour('ca', 'cb', False, '013_contour_plots/figure_DARR_CaCb_bottom.pdf')
# plot_contour('co', 'ca', True,  '013_contour_plots/figure_DARR_COCa_top.pdf')
# plot_contour('co', 'ca', False, '013_contour_plots/figure_DARR_COCa_bottom.pdf')
# plot_contour('co', 'cb', True,  '013_contour_plots/figure_DARR_COCb_top.pdf')
# plot_contour('co', 'cb', False, '013_contour_plots/figure_DARR_COCb_bottom.pdf')
# plot_contour( 'n', 'ca', True,  '013_contour_plots/figure_TEDOR_nCa_nxaxis.pdf')
# plot_contour( 'n', 'cb', True,  '013_contour_plots/figure_TEDOR_nCb_nxaxis.pdf')


if __name__ == "__main__":
  
  main()

