// SWIG Interface to Catalogue

%ignore *::operator[];

%shared_ptr(cbl::catalogue::Object);
%shared_ptr(cbl::catalogue::RandomObject);
%shared_ptr(cbl::catalogue::Halo);
%shared_ptr(cbl::catalogue::Mock);
%shared_ptr(cbl::catalogue::Galaxy);
%shared_ptr(cbl::catalogue::Cluster);
%shared_ptr(cbl::catalogue::Void);
%shared_ptr(cbl::catalogue::Catalogue);

%{
#include "Object.h"
#include "RandomObject.h"
#include "Halo.h"
#include "Mock.h"
#include "Galaxy.h"
#include "Cluster.h"
#include "Catalogue.h"
#include "ChainMesh_Catalogue.h"
#include "Void.h"
%}

%include "Object.h"
%include "RandomObject.h"
%include "Halo.h"
%include "Mock.h"
%include "Galaxy.h"
%include "Cluster.h"
%include "Catalogue.h"
%include "ChainMesh_Catalogue.h"
%include "Void.h"

%template(RandomObjVector) std::vector<cbl::catalogue::RandomObject>;
%template(MockVector) std::vector<cbl::catalogue::Mock>;
%template(HaloVector) std::vector<cbl::catalogue::Halo>;
%template(GalaxyVector) std::vector<cbl::catalogue::Galaxy>;
%template(ClusterVector) std::vector<cbl::catalogue::Cluster>;

%template(VarVector) std::vector<enum cbl::catalogue::Var>;


%extend cbl::catalogue::Catalogue
{  
  %template(add_object) add_object< RandomObject >;
  %template(add_object) add_object< Mock >;
  %template(add_object) add_object< Halo >;
  %template(add_object) add_object< Galaxy >;
  %template(add_object) add_object< Cluster >;
  %template(add_object) add_object< Void >;

  %template(add_objects) add_objects< RandomObject >;
  %template(add_objects) add_objects< Mock >;
  %template(add_objects) add_objects< Halo >;
  %template(add_objects) add_objects< Galaxy >;
  %template(add_objects) add_objects< Cluster >;
  %template(add_objects) add_objects< Void >;
  
  %template(replace_objects) replace_objects< RandomObject >;
  %template(replace_objects) replace_objects< Mock >;
  %template(replace_objects) replace_objects< Halo >;
  %template(replace_objects) replace_objects< Galaxy >;
  %template(replace_objects) replace_objects< Cluster >;
  %template(replace_objects) replace_objects< Void >;

  std::shared_ptr<cbl::catalogue::Object> __getitem__(const size_t i)
    {
      return (*self)[i];
    }
}
