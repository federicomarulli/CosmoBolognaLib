// SWIG Interface to FuncGrid

%module cblFuncGrid

%include "Path.i"
%import "Kernel.i"

%shared_ptr(cbl::glob::FuncGrid);
%shared_ptr(cbl::glob::FuncGrid2D);
%shared_ptr(cbl::glob::FuncGrid_Bspline);

%{
#include "FuncGrid.h"
#include "FuncGrid_Bspline.h"
%}

%include "FuncGrid.h"
%include "FuncGrid_Bspline.h"


%template(FuncGridVector) std::vector<cbl::glob::FuncGrid>;
%template(FuncGrid2DVector) std::vector<cbl::glob::FuncGrid2D>;
%template(FuncGridBsplineVector) std::vector<cbl::glob::FuncGrid_Bspline>;

%template(FuncGridPtrVector) std::vector<std::shared_ptr<cbl::glob::FuncGrid>>;
%template(FuncGrid2DPtrVector) std::vector<std::shared_ptr<cbl::glob::FuncGrid2D>>;
%template(FuncGridBsplinePtrVector) std::vector<std::shared_ptr<cbl::glob::FuncGrid_Bspline>>;
