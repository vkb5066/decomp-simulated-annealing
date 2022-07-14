#!/usr/bin/bash

if [ -z "$1" ]; then
        echo "Enter name of raw output file as a CLA"
        exit 1
fi
echo "Splitting $1 ..."

#Get correct line numbers for splitting from grep
grepOut=(`grep -n "PRINTOUT" $1 | cut -d ':' -f1`)

#Do the actual file splitting
##Enviornment file
echo -n " Chem envs .......... "
tail -n +$(( ${grepOut[0]} + 1)) $1 | head -n $(( ${grepOut[1]} - ${grepOut[0]} - 1)) >$1".env"
echo "written to $1.env"
##Environment decomp file
echo -n " Chem env decomps ... "
tail -n +$(( ${grepOut[2]} + 1)) $1 | head -n $(( ${grepOut[3]} - ${grepOut[2]} - 1)) >$1".dec"
echo "written to $1.dec"
##Site occupancy file
echo -n " Site occupancies ... "
tail -n +$(( ${grepOut[4]} + 1)) $1 | head -n $(( ${grepOut[5]} - ${grepOut[4]} - 1)) >$1".occ"
echo "written to $1.occ"

#Finish up
echo "Done"
