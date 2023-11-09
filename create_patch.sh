#!/bin/bash

# Credits:
#  - `select` bash command based on: https://askubuntu.com/a/1716

###############################################################
## User input & setup
#

printf "\nPlease select the VASP version to patch:\n"
options=("6.4.1" "6.2.1" "5.4.4")
select opt in "${options[@]}"
do
    case $opt in
        "6.4.1")
            printf "you've chosen VASP version 6.4.1\n"
            VDIRNAME="vasp.6.4.1"
            declare -a FILENAMES=(".objects" "main.F" "pot.F" "metagga.F" "nmr.F" "base.F" "wave_struct.F" "reader.F" "xml_writer.F")
            # declare -a FILENAMES=("pot.F")
            break
            ;;
        "6.2.1")
            printf "you've chosen VASP version 6.2.1\n"
            VDIRNAME="vasp.6.2.1"
            declare -a FILENAMES=(".objects" "main.F" "pot.F" "metagga.F" "nmr.F" "base.F" "wave_struct.F" "reader.F" "xml_writer.F")
            # declare -a FILENAMES=("pot.F")
            break
            ;;
        "5.4.4")
            printf "you've chosen VASP version 5.4.4\n"
            VDIRNAME="vasp.5.4.4"
            declare -a FILENAMES=(".objects" "main.F" "pot.F" "metagga.F" "nmr.F" "base.F" "wave.F" "reader.F" "xml_writer.F")
            break
            ;;
        *) echo "invalid option $REPLY";;
    esac
done

PATCHDIRNAME="./${VDIRNAME}/patches"

printf "\nPlease select either the minimal or full patch:\n"
printf "1) minimal patch: applies the Bxc source-free constraint for vasp_ncl with no bells & whistles (no additional INCAR flags or I/O) \n"
printf "2) full patch: the minimal patch, but the SF constraint can be toggled using the LSOURCEFREE flag in the INCAR. In addition, some relevant I/O functionality is included - kindly refer to README file for more information.  \n"
printf "\n"

options=("minimal" "full")
select opt in "${options[@]}"
do
    case $opt in
        "minimal")
            printf "you've chosen to apply the minimal patch\n"
            PATCHDIRNAME="./${VDIRNAME}/patches_minimal"
            declare -a FILENAMES=("pot.F")
            break
            ;;
        "full")
            printf "you've chosen to apply the full patch\n"
            # do nothing
            break
            ;;
        *) echo "invalid option $REPLY";;
    esac
done

printf "\nPlease enter the directory to which you would like to apply the patch:\n"
read -r DIRNAME

###############################################################


###############################################################
###############################################################

# DIRNAME="/global/u1/g/gmoore/conda_code/code/vasp6.2.1/vasp.6.2.1_gpu/src"
# PATCHDIRNAME="./vasp.6.2.1/patches"
# # PATCHDIRNAME="./vasp.6.2.1/patches_buffer"
# declare -a FILENAMES=(".objects" "main.F" "pot.F" "metagga.F" "nmr.F" "base.F" "wave_struct.F" "reader.F" "xml_writer.F")

# DIRNAME="/global/u1/g/gmoore/conda_code/code/vasp.5.4.4/src"
# PATCHDIRNAME="./vasp.5.4.4/patches"
# declare -a FILENAMES=(".objects" "main.F" "pot.F" "metagga.F" "nmr.F" "base.F" "wave.F" "reader.F" "xml_writer.F")

###############################################################

# echo $PATCHDIRNAME
# echo $DIRNAME
# echo ${FILENAMES[0]}

# # exit 1
# echo "GOT HERE"

###############################################################
###############################################################


###############################################################
## Apply/create patches

# # First copy main source file
# cp ${PATCHDIRNAME}/src/source_free_bxc.F ${DIRNAME}

FILETAG_ORIG="orig"
FILETAG_NEW="sf"

# TODO: Currently this file is the same for all versions
# But that won't always be the case, would need an update.
cp vasp.6.2.1/src/source_free_bxc.F "${DIRNAME}/."

for FILENAME in "${FILENAMES[@]}"
do
    
    printf "source code file: ${FILENAME}\n"

    ###############################
    ## Apply patches
    
    # cp ${DIRNAME}/${FILENAME}.orig ${DIRNAME}/${FILENAME}
    
    cp "${DIRNAME}/${FILENAME}" "${DIRNAME}/${FILENAME}.bak"
    patch "${DIRNAME}/${FILENAME}" "${PATCHDIRNAME}/${FILENAME}.patch"
    
    ###############################
    
#     ###############################
#     ## Create patches

#     # diff -Nar -U 0 [...]

#     diff -Nar -U 0 ${DIRNAME}/${FILENAME}.${FILETAG_ORIG} ${DIRNAME}/${FILENAME}.${FILETAG_NEW} > ${PATCHDIRNAME}/${FILENAME}.patch

#     SED_CMD="s+${DIRNAME}/${FILENAME}.${FILETAG_ORIG}+DIR_ORIG/${FILENAME}+g"
#     sed -i ${SED_CMD} ${PATCHDIRNAME}/${FILENAME}.patch

#     SED_CMD="s+${DIRNAME}/${FILENAME}.${FILETAG_NEW}+DIR_NEW/${FILENAME}+g"
#     sed -i ${SED_CMD} ${PATCHDIRNAME}/${FILENAME}.patch

#     sed -i 's/20.* -0.*/2022-08-13 00:00:00.000000000 -0700/g' ${PATCHDIRNAME}/${FILENAME}.patch

#     ###############################
    
    # # to save changes:
    # cp ${DIRNAME}/${FILENAME} ${DIRNAME}/${FILENAME}.sf
    
done

###############################################################
