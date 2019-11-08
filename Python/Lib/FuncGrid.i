// SWIG Interface to FuncGrid

%module cblFuncGrid

%shared_ptr(cbl::glob::FuncGrid);
%shared_ptr(cbl::glob::FuncGrid2D);

%{
#include "FuncGrid.h"
%}

%include "FuncGrid.h"


%template(FuncGridVector) std::vector<cbl::glob::FuncGrid>;
%template(FuncGrid2DVector) std::vector<cbl::glob::FuncGrid2D>;

%template(FuncGridPtrVector) std::vector<std::shared_ptr<cbl::glob::FuncGrid>>;
%template(FuncGrid2DPtrVector) std::vector<std::shared_ptr<cbl::glob::FuncGrid2D>>;
