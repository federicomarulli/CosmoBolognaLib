// SWIG Interface to NumberCounts

%module cblNumberCounts

%import "Path.i"
%import "Kernel.i"
%import "Histogram.i"
%import "Catalogue.i"
%import "Cosmology.i"
%import "Measure.i"

%{
#include "Data.h"
#include "Data1D_collection.h"
#include "CovarianceMatrix.h"
#include "TaperedCovarianceMatrix.h"
#include "CatalogueChainMesh.h"
#include "Chi2.h"
#include "LikelihoodParameters.h"
#include "PosteriorParameters.h"
#include "CombinedPosterior.h"
  
#include "NumberCounts.h"
#include "NumberCounts1D.h"
#include "NumberCounts2D.h"
#include "NumberCounts1D_Redshift.h"
#include "NumberCounts1D_Mass.h"
#include "NumberCounts1D_MassProxy.h"
#include "NumberCounts2D_RedshiftMass.h"
#include "NumberCounts1D_Size.h"
%}

%include "NumberCounts.h"
%include "NumberCounts1D.h"
%include "NumberCounts2D.h"
%include "NumberCounts1D_Redshift.h"
%include "NumberCounts1D_Mass.h"
%include "NumberCounts1D_MassProxy.h"
%include "NumberCounts2D_RedshiftMass.h"
%include "NumberCounts1D_Size.h"
