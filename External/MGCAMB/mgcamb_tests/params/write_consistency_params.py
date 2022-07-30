#
# This file is part of the MGCAMB code
#
# MGCAMB v3.0
#
#   authors:
#       Alex Zucca: azucca@sfu.ca

# This file contains a script that generates a set of CAMB parameters file for testing the GRtransition and the different masses of neutrinos
import pandas as pd
import numpy as np

# this is the list of grtrans and mnu we wish to consider
grtrans_list = [0.0001, 0.001, 0.005, 0.01, 0.05, 0.1]
mnu_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

#grtrans_list = [ 0.001, 0.005, 0.01]
#mnu_list = [0.05, 0.1]


for mnu in mnu_list:
    
    # start by writing the GR file
    file_gr = open('LCDM_GR_mnu_0p'+str(mnu)[2:]+'.ini', 'w')
    omnuh2 = 0.00064/0.06*mnu
    file_gr.write('# omnuh set here \n')
    file_gr.write('omnuh2 = '+str(omnuh2)+'\n')
    file_gr.write('\n')
    file_gr.write('# output_root set here \n')
    file_gr.write('output_root = ./mgcamb_tests/results/results_consistency/LCDM_GR_mnu_0p'+str(mnu)[2:]+'\n')
    file_gr.write('\n')
    file_gr.write('DEFAULT(base_params.ini) \n')
    file_gr.close()
    
    # then for each GR trans generate the MG file
    for grtrans in grtrans_list:
        file =open('LCDM_mg_mnu_0p'+str(mnu)[2:]+'_grt_0p'+str(grtrans)[2:]+'.ini', 'w')
        # just compute the omnuh2 for some mass of neutrinos
        omnuh2 = 0.00064/0.06*mnu
        file.write('# omnuh set here \n')
        file.write('omnuh2 = '+str(omnuh2)+'\n')
        file.write('\n')
        file.write('# output_root set here \n')
        file.write('output_root = ./mgcamb_tests/results/results_consistency/LCDM_mg_mnu_0p'+str(mnu)[2:]+'_grt_0p'+str(grtrans)[2:]+'\n')
        file.write('\n')
        file.write('# MGCAMB flags: \n')
        file.write('MG_flag = 1 \n')
        file.write('pure_MG_flag = 1 \n')
        file.write('GRtrans = '+str(grtrans)+'\n')
        file.write('\n')
        file.write('DEFAULT(base_params.ini) \n')
        file.write('DEFAULT(base_MG_params.ini) \n')
        file.close()





