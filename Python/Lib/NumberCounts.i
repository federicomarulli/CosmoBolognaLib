// SWIG Interface to NumberCounts

%module cblNumberCounts

%{
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
