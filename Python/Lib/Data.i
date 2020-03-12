// SWIG Interface to Data

%module cblData

%shared_ptr(cbl::data::Data);
%shared_ptr(cbl::data::Data1D);
%shared_ptr(cbl::data::Data2D);
%shared_ptr(cbl::data::Data1D_extra);
%shared_ptr(cbl::data::Data2D_extra);
%shared_ptr(cbl::data::Data1D_collection);
%shared_ptr(cbl::data::CovarianceMatrix);
%shared_ptr(cbl::data::TaperedCovarianceMatrix);

%{
#include "Data.h"
#include "Data1D.h"
#include "Data2D.h"
#include "Data1D_extra.h"
#include "Data2D_extra.h"
#include "Data1D_collection.h"
#include "CovarianceMatrix.h"
#include "TaperedCovarianceMatrix.h"
%}

%include "Data.h"
%include "Data1D.h"
%include "Data2D.h"
%include "Data1D_extra.h"
%include "Data2D_extra.h"
%include "Data1D_collection.h"
%include "CovarianceMatrix.h"
%include "TaperedCovarianceMatrix.h"

%template(DataVector) std::vector<cbl::data::Data>;
%template(Data1DVector) std::vector<cbl::data::Data1D>;
%template(Data2DVector) std::vector<cbl::data::Data2D>;
%template(Data1DExtraVector) std::vector<cbl::data::Data1D_extra>;
%template(Data2DExtraVector) std::vector<cbl::data::Data2D_extra>;
%template(CovarianceVector) std::vector<cbl::data::CovarianceMatrix>;
%template(TaperedCovarianceVector) std::vector<cbl::data::TaperedCovarianceMatrix>;

%template(DataPtrVector) std::vector<std::shared_ptr<cbl::data::Data>>;
%template(Data1DPtrVector) std::vector<std::shared_ptr<cbl::data::Data1D>>;
%template(Data2DPtrVector) std::vector<std::shared_ptr<cbl::data::Data2D>>;
%template(Data1DEPtrVector) std::vector<std::shared_ptr<cbl::data::Data1D_extra>>;
%template(Data2DEPtrVector) std::vector<std::shared_ptr<cbl::data::Data2D_extra>>;
%template(CovariancePtrVector) std::vector<std::shared_ptr<cbl::data::CovarianceMatrix>>;
%template(TaperedCovariancePtrVector) std::vector<std::shared_ptr<cbl::data::TaperedCovarianceMatrix>>;
