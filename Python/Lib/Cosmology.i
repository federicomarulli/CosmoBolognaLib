// SWIG Interface to Cosmology

%module cblCosmology

%include "Path.i"
%import "Kernel.i"
%import "Data.i"

%{
#include "Data.h"
#include "Data1D.h"
#include "Data1D_extra.h"
#include "Data1D_collection.h"
#include "TaperedCovarianceMatrix.h"
  
#include "EisensteinHu.h"
#include "CAMB.h"  
#include "Cosmology.h"
#include "LCDM.h"
#include "kLCDM.h"
#include "w0CDM.h"
#include "kw0CDM.h"
#include "w0waCDM.h"
#include "kw0waCDM.h"
#include "LnuCDM.h"
#include "kLnuCDM.h"
#include "w0nuCDM.h"
#include "kw0nuCDM.h"
#include "w0wanuCDM.h"
#include "kw0wanuCDM.h" 
#include "HaloProfile.h"
#include "SuperSampleCovariance.h"
#include "Velocities.h"
#include "3PCF.h"
#include "SizeFunction.h"
#include "BAO.h"
#include "NG.h"
#include "Bias.h"
#include "MassFunction.h"
#include "RSD.h"
#include "MassGrowth.h"
#include "PkXiNonLinear.h"
#include "PkXizSpace.h"
#include "PkXi.h"
#include "Sigma.h"
%}

%include <std_shared_ptr.i>
%shared_ptr(cbl::cosmology::Cosmology)
%shared_ptr(cbl::cosmology::LCDM)
%shared_ptr(cbl::cosmology::kLCDM)
%shared_ptr(cbl::cosmology::w0CDM)
%shared_ptr(cbl::cosmology::kw0CDM)
%shared_ptr(cbl::cosmology::w0waCDM)
%shared_ptr(cbl::cosmology::kw0waCDM)
%shared_ptr(cbl::cosmology::LnuCDM)
%shared_ptr(cbl::cosmology::kLnuCDM)
%shared_ptr(cbl::cosmology::w0nuCDM)
%shared_ptr(cbl::cosmology::kw0nuCDM)
%shared_ptr(cbl::cosmology::w0wanuCDM)
%shared_ptr(cbl::cosmology::kw0wanuCDM) 

%include "Data.h"
%include "Data1D.h"
%include "Data1D_extra.h"
%include "Data1D_collection.h"
%include "TaperedCovarianceMatrix.h"
  
%include "EisensteinHu.h"
%include "CAMB.h"  
%include "Cosmology.h"
%include "LCDM.h"
%include "kLCDM.h"
%include "w0CDM.h"
%include "kw0CDM.h"
%include "w0waCDM.h"
%include "kw0waCDM.h"
%include "LnuCDM.h"
%include "kLnuCDM.h"
%include "w0nuCDM.h"
%include "kw0nuCDM.h"
%include "w0wanuCDM.h"
%include "kw0wanuCDM.h" 
%include "HaloProfile.h"
%include "SuperSampleCovariance.h"
%include "Velocities.h"
%include "3PCF.h"
%include "SizeFunction.h"
%include "BAO.h"
%include "NG.h"
%include "Bias.h"
%include "MassFunction.h"
%include "RSD.h"
%include "MassGrowth.h"
%include "PkXiNonLinear.h"
%include "PkXizSpace.h"
%include "PkXi.h"
%include "Sigma.h"
