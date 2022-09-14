

# DIRNAME="/global/u1/g/gmoore/conda_code/code/vasp6.2.1/vasp.6.2.1_gpu/src"
# PATCHDIRNAME="./vasp.6.2.1/patches"
# # PATCHDIRNAME="./vasp.6.2.1/patches_buffer"
# declare -a FILENAMES=(".objects" "main.F" "pot.F" "metagga.F" "nmr.F" "base.F" "wave_struct.F" "reader.F" "xml_writer.F")


DIRNAME="/global/u1/g/gmoore/conda_code/code/vasp.5.4.4/src"
PATCHDIRNAME="vasp.5.4.4/patches"
declare -a FILENAMES=(".objects" "main.F" "pot.F" "metagga.F" "nmr.F" "base.F" "wave.F" "reader.F" "xml_writer.F")


# # First copy main source file
# cp vasp.6.2.1/src/source_free_bxc.F ${DIRNAME}

FILETAG_ORIG="orig"
FILETAG_NEW="sf"

for FILENAME in "${FILENAMES[@]}"
do
   
#     ###############################
#     ## Create patches

#     # diff -Naru [...]

#     diff -Nar -U 0 ${DIRNAME}/${FILENAME}.${FILETAG_ORIG} ${DIRNAME}/${FILENAME}.${FILETAG_NEW} > ${PATCHDIRNAME}/${FILENAME}.patch

#     SED_CMD="s+"${DIRNAME}"/"${FILENAME}"."${FILETAG_ORIG}"+DIR_ORIG/"${FILENAME}"+g"
#     sed -i ${SED_CMD} ${PATCHDIRNAME}/${FILENAME}.patch

#     SED_CMD="s+"${DIRNAME}"/"${FILENAME}"."${FILETAG_NEW}"+DIR_NEW/"${FILENAME}"+g"
#     sed -i ${SED_CMD} ${PATCHDIRNAME}/${FILENAME}.patch

#     sed -i 's/20.* -0.*/2022-08-13 00:00:00.000000000 -0700/g' ${PATCHDIRNAME}/${FILENAME}.patch

#     ###############################


    ###############################
    ## Apply patches
    
    cp ${DIRNAME}/${FILENAME}.orig ${DIRNAME}/${FILENAME}
    
    patch ${DIRNAME}/${FILENAME} ${PATCHDIRNAME}/${FILENAME}.patch
    
    ###############################
    
    # # to save changes:
    # cp ${DIRNAME}/${FILENAME} ${DIRNAME}/${FILENAME}.sf
    
done

