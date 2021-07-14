#----------------------------------------------------------------------------------------
#
# This file is part of MGCAMB.
#
#----------------------------------------------------------------------------------------

#
# This file contains the make targets for the test suite
#

# structure of the folder:

TEST_DIR=$(abspath .)
SCRIPT_DIR=$(TEST_DIR)/test_scripts
RESULTS_DIR=$(TEST_DIR)/results
PARAMETER_DIR=$(TEST_DIR)/parameters
PYTHON_DIR=$(TEST_DIR)/python
LEGACY_DIR=$(TEST_DIR)/eftcamb_legacy

# make general targets:

all: consistency consistency_plots

test: consistency

default: test

# results creation target:

consistency:
	@bash $(SCRIPT_DIR)/test_consistency.sh


# plot target:

consistency_plots:
	@bash $(SCRIPT_DIR)/consistency_plot.sh
	
# clean target:

clean:
	@rm -rf $(RESULTS_DIR)/spectra_results/*
	@rm -rf $(RESULTS_DIR)/spectra_plots/*
	@rm -rf $(RESULTS_DIR)/profile/*
	@rm -rf $(RESULTS_DIR)/legacy_spectra_plot/*
	@rm -rf $(RESULTS_DIR)/benchmark_results/*
