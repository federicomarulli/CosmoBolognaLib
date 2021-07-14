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
printf "MGCAMB TEST: producing consisteny plots\n"
printf "********************************************\n${Color_Off}"

# go to the python directory:
cd $TEST_PYTHON_DIR

SUCCESS=true

printf "  Making plots "
python ./consistency_plots.py  &> $RESULTS_LOGS/consistency_plots.log
python ./mpk_plot.py  &> $RESULTS_LOGS/mpk_plot.log

# check if run succeded:
if [ $? -eq 0 ]; then
    printf "${BGreen} OK\n${Color_Off}"
else
    printf "${BRed} FAIL\n${Color_Off}"
    SUCCESS=false
fi

printf "\n"

if [ "$SUCCESS" = true ]; then
	    printf "${Green}********************************************\n"
		printf "MGCAMB TEST: done consistency plots\n"
		printf "********************************************\n${Color_Off}"
else
	    printf "${Yellow}********************************************\n"
		printf "MGCAMB TEST: done consistency plots\n"
		printf "WARNING: plotting script failed.\n"
		printf "********************************************\n${Color_Off}"
fi


# return status:

if [ "$SUCCESS" = true ]; then
	    exit 0
else
	    exit 1
fi
