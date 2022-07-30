#!/usr/bin/env python

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

import sys
import os

# import CosmoBolognaLib functions #
import CosmoBolognaLib as cbl


##########################################################
# just remember that parameter files do have a purpose

filename = "../input/parameter_file.ini"
print(" Loading parameters from", filename)
param = cbl.ReadParameters(filename)
  
##########################################################
# Coordinates type selection

if param.findBool('comovingCoordinates') :
  coordinates = cbl.CoordinateType__comoving_
else :
  coordinates = cbl.CoordinateType__observed_

##########################################################
# load the input void catalogue

cast = []
clmn = []
attrNames = ['X_coord', 'Y_coord', 'Z_coord', 'Radius', 'centralDensity', 'densityContrast']
attrAv = [cbl.Var__X_, cbl.Var__Y_, cbl.Var__Z_, cbl.Var__Radius_, cbl.Var__CentralDensity_, cbl.Var__DensityContrast_]
for ii in range(len(attrNames)) :
  if param.findBool(attrNames[ii]) :
    cast.append(attrAv[ii])
    clmn.append(param.findInt(attrNames[ii]+'_clmn'))
clmn, cast = (list(x) for x in zip(*sorted(zip(clmn, cast))))  # orders clmn and cast according to column order 

attr = cbl.VarCast(cast)

vdcat = cbl.Catalogue (cbl.ObjectType__Void_,
                       coordinates,
                       attr,
                       clmn,
                       [param.findString('inputVoidCatalogue')],
                       param.findInt('vd_comments'))

if (param.findDouble('boxside') < 0.) :
  boxside = abs(vdcat.Max(cbl.Var__X_) - vdcat.Min(cbl.Var__X_))
else :
  boxside = param.findDouble('boxside')

      
##########################################################
# load the input tracers catalogue

if param.findBool('Gadget') :
  if not param.findBool('comovingCoordinates') :
    ErrorCBL('Observed coordinates not available for Gadget snapshot.')
  else :
    trcat = cbl.Catalogue (cbl.ObjectType__Halo_,
                           param.findString('inputTracersFile'),
                           param.findBool('swapEndianism'),
                           param.findDouble('fact'),
                           True,
                           param.findDouble('nSub'))
else :
  if param.findBool('comovingCoordinates') :
    
    tr_cast = []
    tr_clmn = []
    trAttrNames = ['X_coord_tr', 'Y_coord_tr', 'Z_coord_tr', 'Mass']
    trAttrAv = [cbl.Var__X_, cbl.Var__Y_, cbl.Var__Z_, cbl.Var__Mass_]
    for ii in range(len(trAttrNames)) :
      if param.findBool(trAttrNames[ii]) :
        tr_cast.append(trAttrAv[ii])
        tr_clmn.append(param.findInt(trAttrNames[ii]+'_clmn'))        
    tr_clmn, tr_cast = (list(x) for x in zip(*sorted(zip(tr_clmn, tr_cast))))  # orders clmn and cast according to column order 
    tr_attr = cbl.VarCast(tr_cast)
    
    temp = cbl.Catalogue (cbl.ObjectType__Halo_,
                          coordinates,
                          tr_attr,
                          tr_clmn,
                          [param.findString('inputTracersFile')],
                          param.findInt('tr_comments'),
                          param.findDouble('nSub'),
                          param.findDouble('fact'))

    if not param.findBool('Mass') :
      
      trcat = temp
      temp = None
      print(" Finished reading input tracers catalogue.")
      
    else :
      
      print("Finished reading input tracers catalogue, now applying mass scale factor and/or cut-off ... ")

      # scale factor
      if (param.findDouble('Munit') > 0.) :
        for ii in range(temp.nObjects()) :
          mass = temp.mass(ii)*param.findDouble('Munit')
          temp.set_var(ii, cbl.Var__Mass_, mass)

      # mass cut-off
      if (param.findDouble('Mmin') > 0.) :
        trcat = cbl.Catalogue ()
        trcat = temp.sub_catalogue(cbl.Var__Mass_, param.findDouble('Mmin'), temp.Max(cbl.Var__Mass_), False)
      else :
        trcat = temp
        temp = None

      print("\t ... done!")

  # observed coordinates
  else :
    print("Observed coordinates not supported with this configuration...")
    exit(1)

  
  
##########################################################
# Finally building the cleaned catalogue

# sets the radius if not read from file:
if not param.findBool('Radius') :
  delta_r = param.findVectorDouble('delta_r')
  radii = [delta_r[1] for ii in range(vdcat.nObjects())]
  vdcat.set_var(cbl.Var__Radius_, radii)

# Generate chain-mesh of the input tracers catalogue              
ChM = cbl.ChainMesh3D (2.*trcat.mps(),
                       trcat.var(cbl.Var__X_),
                       trcat.var(cbl.Var__Y_),
                       trcat.var(cbl.Var__Z_),
                       vdcat.Max(cbl.Var__Radius_))

# sets the central density if not read from file:
#if not param.findBool('centralDensity') : 
#  vdcat.compute_centralDensity(trcat,ChM)

# sets the density contrast if not read from file:
#if not param.findBool('densityContrast') :
#  vdcat.compute_densityContrast(trcat,ChM,param.findDouble('ratio'))

# overlap-check criterion choice:
ol_crit = cbl.Var__DensityContrast_ if param.findInt('ol_crit') == 1 else cbl.Var__CentralDensity_

# choose the (positive) value of the threshold
threshold = 1. + param.findDouble('deltav_NL')

# build the catalogue:
vdcat.clean_void_catalogue(	param.findBool('initial_radius'),
                               param.findVectorDouble('delta_r'),
                               threshold,
                               param.findBool('rescale'),
                               trcat,
                               ChM,
                               param.findDouble('ratio'),
                               param.findBool('overlap'),
                               ol_crit)
print('\n')

# write the obtained catalogue in a file

clmnsToPrint = cbl.VarCast(attrAv)

if not os.path.exists(param.findString('outputDir')):
    os.makedirs(param.findString('outputDir'))

vdcat.write_data(param.findString('outputDir')+param.findString('VoidCatalogueFile'),
                         clmnsToPrint)

#############################################################################################################
##################################### .. and that's all folks!! #############################################
#############################################################################################################
