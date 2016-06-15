// SWIG Interface to Catalogue

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

%template(RandomObjVec) vector<cosmobl::catalogue::RandomObject>;
%template(MockVec) vector<cosmobl::catalogue::Mock>;
%template(HaloVec) vector<cosmobl::catalogue::Halo>;
%template(GalaxyVec) vector<cosmobl::catalogue::Galaxy>;
%template(ClusterVec) vector<cosmobl::catalogue::Cluster>;
%template(ObjectPtr) shared_ptr<cosmobl::catalogue::Object>;
%template(VoidPtr) shared_ptr<cosmobl::catalogue::Void>;

%extend cosmobl::catalogue::Catalogue
{  
  %template(Catalogue) Catalogue< RandomObject >;
  %template(Catalogue) Catalogue< Mock >;
  %template(Catalogue) Catalogue< Halo >;
  %template(Catalogue) Catalogue< Galaxy >;
  %template(Catalogue) Catalogue< Cluster >;
  //%template(Catalogue) Catalogue< Void >;
 
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
  
  %template(remove_objects) remove_objects< RandomObject >;
  %template(remove_objects) remove_objects< Mock >;
  %template(remove_objects) remove_objects< Halo >;
  %template(remove_objects) remove_objects< Galaxy >;
  %template(remove_objects) remove_objects< Cluster >;
  %template(remove_objects) remove_objects< Void >;
}

