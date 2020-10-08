// SWIG Interface to CosmoBolognaLib

%module(directors="1") CosmoBolognaLib

//%import "../Doc/documentation.i"

%include Lib/Kernel.i
%include Lib/Wrappers.i
%include Lib/FuncGrid.i
%include Lib/FFT.i
%include Lib/Random.i
%include Lib/Func.i
%include Lib/Data.i
%include Lib/Field.i
%include Lib/Histogram.i
%include Lib/Distribution.i
%include Lib/Stat.i
%include Lib/Cosmology.i
%include Lib/ChainMesh.i
%include Lib/Catalogue.i
%include Lib/LogNormal.i
%include Lib/Measure.i
%include Lib/NumberCounts.i
%include Lib/TwoPointCorrelation.i
%include Lib/ThreePointCorrelation.i
%include Lib/CovarianceMatrix.i
%include Lib/Modelling.i
%include Lib/Modelling_Cosmology.i
%include Lib/Modelling_NumberCounts.i
%include Lib/Modelling_TwoPointCorrelation.i
%include Lib/Modelling_ThreePointCorrelation.i
%include Lib/GlobalFunc.i
%include Lib/ReadParameters.i
