// SWIG Interface to Catalogue

%module cblCatalogue

%{
#include "../../Headers/Objects/Object.h"
#include "../../Headers/Objects/GenericObject.h"
#include "../../Headers/Objects/Halo.h"
#include "../../Headers/Objects/Mock.h"
#include "../../Headers/Objects/Galaxy.h"
#include "../../Headers/Objects/Cluster.h"
#include "../../Headers/Lib/Catalogue.h"
#include "../../Headers/Lib/ChainMesh_Catalogue.h"
%}

%include "../../Headers/Objects/Object.h"
%include "../../Headers/Objects/GenericObject.h"
%include "../../Headers/Objects/Halo.h"
%include "../../Headers/Objects/Mock.h"
%include "../../Headers/Objects/Galaxy.h"
%include "../../Headers/Objects/Cluster.h"
%include "../../Headers/Lib/Catalogue.h"
%include "../../Headers/Lib/ChainMesh_Catalogue.h"
