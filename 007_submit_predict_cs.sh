#!/bin/bash

### run as a bash script
### or submit as an sbatch job

#SBATCH --job-name=ppm
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --time=00-04:00:00
#SBATCH --output=single.%j.out
#SBATCH --error=single.%j.err

PPM_executable="/home/user/foo/PPM/ppm_linux"
# script assumes the input folders are named: 005_X_fm2_fixbb
#   and that the output is into folders are named: 008_X_fm2_cs

for input_folder in 005_*_fm2_fixbb ; do
    ## figure out in which folder the output will go
    remove_fileend=${input_folder%_fm2_fixbb}
    amino_acid=${remove_fileend#005_}
    output_folder="008_${amino_acid}_fm2_cs"

    ## run PPM for each PDB file
    ## output the .cs files into the right output folder
    for pdb_path in $( ls -1 $input_folder/*.pdb ); do
        pdb_file=${pdb_path#*/}
        filename_core=${pdb_file%.pdb}
#       echo ${input_folder}/${pdb_file} " -> " ${output_folder}/${filename_core}.cs
        $PPM_executable -para old -pdb ${input_folder}/${pdb_file} -pre ${output_folder}/${filename_core}.cs
    done

    ## cleanup
    rm -f bb_details.dat
    rm -f bb_predict.dat
    rm -f cs_rmsd.dat
    rm -f proton_details.dat
    rm -f proton_predict.dat
done
