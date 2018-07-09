// SWIG Interface to Data

%module cblData

%shared_ptr(cbl::data::Data);
%shared_ptr(cbl::data::Data1D);
%shared_ptr(cbl::data::Data2D);
%shared_ptr(cbl::data::Data1D_extra);
%shared_ptr(cbl::data::Data2D_extra);
%shared_ptr(cbl::data::Data1D_collection);

%{
#include "Data.h"
#include "Data1D.h"
#include "Data2D.h"
#include "Data1D_extra.h"
#include "Data2D_extra.h"
#include "Data1D_collection.h"
%}

%include "Data.h"
%include "Data1D.h"
%include "Data2D.h"
%include "Data1D_extra.h"
%include "Data2D_extra.h"
%include "Data1D_collection.h"
