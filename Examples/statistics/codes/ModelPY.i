%module modelpy

%include "stl.i"
%include "stdint.i"
%include "std_string.i"
%include "std_vector.i"
%include <std_shared_ptr.i>

%shared_ptr(cbl::statistics::Model1D);
%shared_ptr(cbl::cosmology::Cosmology);

%import (module="CosmoBolognaLib") "Model1D.h"
%import (module="CosmoBolognaLib") "Cosmology.h"

%{
#include "ModelPY.h"
%}

%include "ModelPY.h"
