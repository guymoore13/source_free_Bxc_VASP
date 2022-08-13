
DIFFPATCHDIRNAME="./vasp.5.4.4/diff_patches_6.2.1"
PATCHDIRNAME="./vasp.5.4.4/patches"

declare -a FILENAMES=("pot.F")
# declare -a FILENAMES=(".objects" "main.F" "pot.F" "metagga.F" "nmr.F" "base.F" "wave_struct.F" "reader.F" "xml_writer.F")

for FILENAME in "${FILENAMES[@]}"
do
   
   patch -p0 -i ${DIFFPATCHDIRNAME}/${FILENAME}.diff.patch -o ${PATCHDIRNAME}/${FILENAME}.patch
   
done
