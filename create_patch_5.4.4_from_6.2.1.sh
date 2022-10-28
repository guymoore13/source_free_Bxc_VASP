#!/bin/bash

DIFFPATCHDIRNAME="./vasp.5.4.4/diff_patches_6.2.1"
PATCHDIRNAME_A="./vasp.5.4.4/patches"
PATCHDIRNAME_B="./vasp.6.2.1/patches"

declare -a FILENAMES_A=(".objects" "main.F" "pot.F" "metagga.F" "nmr.F" "base.F" "wave.F" "reader.F" "xml_writer.F")
declare -a FILENAMES_B=(".objects" "main.F" "pot.F" "metagga.F" "nmr.F" "base.F" "wave_struct.F" "reader.F" "xml_writer.F")

for i in "${!FILENAMES_A[@]}"
do
    
    FILENAME_A="${FILENAMES_A[i]}"
    FILENAME_B="${FILENAMES_B[i]}"
    
    
#     ###############################
#     ## Create diff of patches
    
#     diff -Nar -U 0 ${PATCHDIRNAME_B}/${FILENAME_B}.patch ${PATCHDIRNAME_A}/${FILENAME_A}.patch > ${DIFFPATCHDIRNAME}/${FILENAME_A}.diff.patch
#     sed -i 's/20.* -0700/2022-08-13 00:00:00.000000000 -0700/g' ${DIFFPATCHDIRNAME}/${FILENAME_A}.diff.patch
    
#     ###############################
    
    
    ###############################
    ## Apply diff of patches
    
    patch -p0 -i ${DIFFPATCHDIRNAME}/${FILENAME_A}.diff.patch -o ${PATCHDIRNAME_A}/${FILENAME_A}.patch ${PATCHDIRNAME_B}/${FILENAME_B}.patch
    
    ###############################
    
done
