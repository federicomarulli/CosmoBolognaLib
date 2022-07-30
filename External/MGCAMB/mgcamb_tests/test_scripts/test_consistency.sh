#----------------------------------------------------------------------------------------
#
# This file is part of MGCAMB.
#
#----------------------------------------------------------------------------------------
#
# Script that runs the GR consistency tests of MGCAMB.
#
# Developed by: Alex Zucca (azucca@sfu.ca) for the MGCAMB code
#
#      based on the test suite of the EFTCAMB code developed by Marco Raveri
#

#!/bin/bash

# get the path of the script:
SCRIPT_PATH="`dirname \"$0\"`"                  # relative
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$SCRIPT_PATH" ] ; then
  exit 1
fi

# import things:
source $SCRIPT_PATH/paths.sh
source $SCRIPT_PATH/colors.sh

# start the script:
printf "${Green}********************************************\n"
printf "MGCAMB TEST: running GR consistency checks\n"
printf "********************************************\n${Color_Off}"

# go to the Paramaters directory and generate the param files
cd $TEST_PARAMS_DIR

printf "\nGenerating params files: "
python write_consistency_params.py

# check if compilation succeded:
if [ $? -eq 0 ]; then
    printf "${BGreen} OK\n${Color_Off}"
else
    printf "${BRed} FAIL\n${Color_Off}"
    exit 1
fi


# go to the CAMB directory:
cd $CAMB_DIR

printf "\nCompiling MGCAMB: "

make clean &> /dev/null
make camb &> $RESULTS_LOGS/GR_consistency_log.txt

# check if compilation succeded:
if [ $? -eq 0 ]; then
    printf "${BGreen} OK\n${Color_Off}"
else
    printf "${BRed} FAIL\n${Color_Off}"
    exit 1
fi

printf "\n"

# define flag to check wether all went correct:
SUCCESS=true

# run all the parameters:

for i in $TEST_PARAMS_DIR/*.ini;
    do

    filename=$(basename "$i")
    extension="${filename##*.}"
    filename="${filename%.*}"

    printf "  Doing %s: " "$filename"

    ./camb $i &> $RESULTS_LOGS/$filename.log

    # check if run succeded:
    if [ $? -eq 0 ]; then
        printf "${BGreen} OK\n${Color_Off}"
    else
        printf "${BRed} FAIL\n${Color_Off}"
        SUCCESS=false
    fi

done;

printf "\n"


# cleaning up:
make clean        &> /dev/null
rm   camb  &> /dev/null

# feedback for the run:

printf "Results in: %s \n" "$RESULTS_CONSISTENCY_DIR"
printf "\n"

if [ "$SUCCESS" = true ]; then
	    printf "${Green}********************************************\n"
		printf "MGCAMB TEST: done consistency runs\n"
		printf "********************************************\n${Color_Off}"
else
	    printf "${Yellow}********************************************\n"
		printf "MGCAMB TEST: done consistency runs\n"
		printf "WARNING: some parameters failed.\n"
		printf "********************************************\n${Color_Off}"
fi

# produce the global benchmark report:


# return status:

if [ "$SUCCESS" = true ]; then
	    exit 0
else
	    exit 1
fi
