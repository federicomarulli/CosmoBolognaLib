#!/usr/bin/env python

import sys

# import CosmoBolognaLib functions #
import CosmoBolognaLib as cbl                                 
from CosmoBolognaLib import EnumTypes as et
from CosmoBolognaLib import ErrorCBL

# import timer
import time

##########################################################
# just remember that parameter files do have a purpose

if (len(sys.argv) == 1):
  print "Usage: /path/to/cleaVoidCatalogue.py /path/to/parameter_file.ini"
  exit(-1)

if (len(sys.argv) > 1):
  filename = sys.argv[1]
  print " Loading parameters from", filename
  param = cbl.ReadParameters(filename)

else:
  print " Using default parameters"

  
##########################################################
# Coordinates type selection

if param.findBool('comovingCoordinates') :
  coordinates = et._comovingCoordinates_
else :
  coordinates = et._observedCoordinates_

##########################################################
# define a cosmological model, using parameters from file #

cosm = cbl.Cosmology(param.findDouble('OmM'),
                     param.findDouble('Omb'),
                     param.findDouble('Omn'),
                     param.findDouble('massless'),
                     param.findInt('massive'),
                     param.findDouble('OmL'),
                     param.findDouble('Omr'),
                     param.findDouble('hh'),
                     param.findDouble('As'),
                     param.findDouble('pivot'),
                     param.findDouble('ns'),
                     param.findDouble('w0'),
                     param.findDouble('wa'),
                     param.findDouble('fNL'),
                     param.findInt('type_NG'),
                     param.findDouble('tau'),
                     param.findString('model'),
                     param.findBool('unit'))

# set sigma_8 value from file
cosm.set_sigma8(param.findDouble('sigma8'))

##########################################################
# load the input void catalogue

cast = []
clmn = []
attrNames = ['X_coord', 'Y_coord', 'Z_coord', 'Radius', 'centralDensity', 'densityContrast']
attrAv = [et._X_, et._Y_, et._Z_, et._Radius_, et._CentralDensity_, et._DensityContrast_]
for ii in range(len(attrNames)) :
  if param.findBool(attrNames[ii]) :
    cast.append(attrAv[ii])
    clmn.append(param.findInt(attrNames[ii]+'_clmn'))
clmn, cast = (list(x) for x in zip(*sorted(zip(clmn, cast))))  # orders clmn and cast according to column order 
attr = cbl.VarVec()
cbl.VarCast(attr, cast)

vdcat = cbl.Catalogue (et._Void_,
                       coordinates,
                       attr,
                       clmn,
                       [param.findString('inputVoidCatalogue')],
                       param.findInt('vd_comments'))

if (param.findDouble('boxside') < 0.) :
  boxside = abs(vdcat.Max(et._X_) - vdcat.Min(et._X_))
else :
  boxside = param.findDouble('boxside')
vdcat.compute_catalogueProperties(param.findDouble('boxside'))

##########################################################
# load the input tracers catalogue

if param.findBool('Gadget') :
  if not param.findBool('comovingCoordinates') :
    ErrorCBL('Observed coordinates not available for Gadget snapshot.')
  else :
    trcat = cbl.Catalogue (et._Halo_,
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
    trAttrAv = [et._X_, et._Y_, et._Z_, et._Mass_]
    for ii in range(len(trAttrNames)) :
      if param.findBool(trAttrNames[ii]) :
        tr_cast.append(trAttrAv[ii])
        tr_clmn.append(param.findInt(trAttrNames[ii]+'_clmn'))
        
    tr_clmn, tr_cast = (list(x) for x in zip(*sorted(zip(tr_clmn, tr_cast))))  # orders clmn and cast according to column order 
    tr_attr = cbl.VarVec()
    cbl.VarCast(tr_attr, tr_cast)
    
    temp = cbl.Catalogue (et._Halo_,
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
      print "Finished reading input tracers catalogue."
      
    else :
      
      print "Finished reading input tracers catalogue, now sorting and applying mass scale factor and/or cut-off ... "

      # scale factor
      for ii in range(temp.nObjects()) :
        mass = temp.mass(ii)*param.findDouble('Munit')
        temp.set_var(ii, et._Mass_, mass)

      # sorting
      temp.sort(et._Mass_, False)

      # mass cut-off
      if (param.findDouble('Mmin') > 0.) :
        trcat = cbl.Catalogue ()
        ii = 0
        while (temp.mass(ii) > param.findDouble('Mmin')) :
          trcat.add_object(temp.catalogue_object(ii))
          ii = ii + 1
      else :
        trcat = temp
      
      temp = None
      print "\t ... done!"

  # observed coordinates
  else :
    print "Observed coordinates not supported yet..."
    exit(1)

trcat.compute_catalogueProperties(param.findDouble('boxside'))
  
##########################################################
# Generate chain-mesh of the input tracers catalogue              

ChM = cbl.ChainMesh3D (2.*trcat.mps(),
                       trcat.var(et._X_),
                       trcat.var(et._Y_),
                       trcat.var(et._Z_),
                       vdcat.Max(et._Radius_))
  
##########################################################
# Finally building the cleaned catalogue

# sets the radius if not read from file:
if not param.findBool('Radius') :
  limit = param.findVectorDouble('delta_r')
  radii = [delta_r[1] for ii in range(vdcat.nObjects())]
  vdcat.set_var(et._Radius_, radii)

# sets the central density if not read from file:
if not param.findBool('centralDensity') : 
  vdcat.compute_centralDensity(trcat,
                               ChM,
                               trcat.numdensity(),
                               param.findDouble('ratio'))

# sets the density contrast if not read from file:
if not param.findBool('densityContrast') :
  vdcat.compute_densityContrast(trcat,
                                ChM,
                                param.findDouble('ratio'))

# overlap-check criterion choice:
ol_crit = et._DensityContrast_ if param.findInt('ol_crit') == 1 else et._CentralDensity_

# build the catalogue:
print '\n'
tw0 = time.time()
tc0 = time.clock()
vdcat_cleaned = cbl.Catalogue (vdcat,
                               [param.findBool('clean1'),
                                param.findBool('clean2'),
                                param.findBool('clean3')],
                               param.findVectorDouble('delta_r'),
                               param.findDouble('threshold'),
                               param.findDouble('relevance'),
                               param.findBool('rescale'),
                               trcat,
                               ChM,
                               param.findDouble('ratio'),
                               param.findBool('overlap'),
                               ol_crit)
print 'Cleaning the catalogue took: ', time.clock()-tc0, ' sec'
print 'Wall time: ', time.time()-tw0, ' sec'
print '\n'

# write the obtained catalogue in a file

clmnsToPrint = cbl.VarVec()
cbl.VarCast(clmnsToPrint, attrAv)
vdcat_cleaned.write_data(param.findString('outputFile'),
                         clmnsToPrint)

#############################################################################################################
##################################### .. and that's all folks!! #############################################
#############################################################################################################
