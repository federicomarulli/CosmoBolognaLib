#----------------------------------------------------------------------------------------
#
# This file is part of EFTCAMB.
#
# Copyright (C) 2013-2017 by the EFTCAMB authors
#
# The EFTCAMB code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file eftcamb/LICENSE at
# the top level of the EFTCAMB distribution.
#
#----------------------------------------------------------------------------------------

#
# This file contains the paths of the EFTCAMB test suite.
#
# Developed by: Marco Raveri (mraveri@sissa.it) for the EFTCAMB code
#

# get the path of the script:
SCRIPT_PATH="`dirname \"$0\"`"                  # relative
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$SCRIPT_PATH" ] ; then
  exit 1
fi

# get the path of the things in the test suite:

# test directory:
TEST_DIR=$SCRIPT_PATH/..

# camb directory:
CAMB_DIR=$TEST_DIR/..

# test directory:
TEST_PARAMS_DIR=$TEST_DIR/params
TEST_PYTHON_DIR=$TEST_DIR/python
TEST_RESULTS_DIR=$TEST_DIR/results

# results directory:
RESULTS_CONSISTENCY_DIR=$TEST_RESULTS_DIR/results_consistency
RESULTS_CONSISTENCY_PLOTS=$TEST_RESULTS_DIR/consistency_plots
RESULTS_LOGS=$TEST_RESULTS_DIR/logs

