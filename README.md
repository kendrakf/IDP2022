_Code from the article:_
# Peak shapes under DNP-conditions report on conformation of the ensemble
_authors:_ Jaka Kragelj, Rania Dumarieh, Yiling Xiao, Kendra K. Frederick

Summary of contents: 
| File                                | Description                                                                       |
| :---------------------------------- | :-------------------------------------------------------------------------------- |
| 005_PDBs.tar.gz                     | Folder containing PDB files of the Ala(25)-X-Ala(25) ensembles                    |
| 006_make_cs_folders.sh              | Script that creates the 008_PPM_pred folders (in this case they already exist).   |
| 007_submit_predict_cs.sh            | Script for running/submitting chemical shift calculations with PPM_One.           |
| 008_PPM_pred.tar.gz                 | Folder containing chemical shift predictions for the Ala(25)-X-Ala(25) ensembles. |
| 009_make_joined_rama_cs_folders.sh  | Script that creates the 011_PPM_PDB_join folders (here, they already exist).      |
| 010_run_all_join_dihedral_cs.py     | Calculate dihedral angles from PDBs and append the chemical shift information.    |
| 011_PDB_PPM_join.tar.gz             | Folder containing files with phi/psi and chemical shift information for each PDB. |
| 012_run_all_contour_plots.py        | Script for generating the contour plots. |
| 013_contour_plots                   | Folder with the contour plots imitating peak shapes. |
| 014_run_all_1D_histograms.py        | Script for generating the 1D histograms with bars colored according to region.    |
| 015_1D_histogram_plots              | Folder with 1D histograms with bars colored according to regions.                 |
| 016_run_add_gaussian_linewidth.py   | Script for 1D peak shapes (sum of gaussians for every chemical shift prediction). |
| 017_gaussian_simulations            | Folder with 1D peak shape predictions.                                            |

