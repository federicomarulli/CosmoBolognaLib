// SWIG Interface to Modelling_NumberCounts

%module cblModelling_NumberCounts

%include "Path.i"
%include "Kernel.i"
%import "Cosmology.i"
%include "Modelling.i"

%{
#include "Data.h"
#include "Data1D.h"
#include "Data1D_collection.h"
#include "Data1D_extra.h"
#include "CovarianceMatrix.h"
#include "TaperedCovarianceMatrix.h"
#include "Chi2.h"
#include "LikelihoodParameters.h"
#include "PosteriorParameters.h"
#include "CombinedPosterior.h"

#include "ModelFunction_NumberCounts.h"
#include "ModelFunction_NumberCounts1D_Redshift.h"
#include "ModelFunction_NumberCounts1D_Mass.h"
#include "ModelFunction_NumberCounts1D_MassProxy.h"
#include "ModelFunction_NumberCounts2D_RedshiftMass.h"
#include "Modelling_NumberCounts.h"
#include "Modelling_NumberCounts1D.h"
#include "Modelling_NumberCounts2D.h"
#include "Modelling_NumberCounts1D_Redshift.h"
#include "Modelling_NumberCounts1D_Mass.h"
#include "Modelling_NumberCounts1D_MassProxy.h"
#include "Modelling_NumberCounts2D_RedshiftMass.h"
%}

%include "ModelFunction_NumberCounts.h"
%include "ModelFunction_NumberCounts1D_Redshift.h"
%include "ModelFunction_NumberCounts1D_Mass.h"
%include "ModelFunction_NumberCounts1D_MassProxy.h"
%include "ModelFunction_NumberCounts2D_RedshiftMass.h"
%include "Modelling_NumberCounts.h"
%include "Modelling_NumberCounts1D.h"
%include "Modelling_NumberCounts2D.h"
%include "Modelling_NumberCounts1D_Redshift.h"
%include "Modelling_NumberCounts1D_Mass.h"
%include "Modelling_NumberCounts1D_MassProxy.h"
%include "Modelling_NumberCounts2D_RedshiftMass.h"
