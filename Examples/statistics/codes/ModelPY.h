// ====================================================================================
// Example code: write a C++ function and create a class inheriting from Model1D object
// ====================================================================================

#include "Cosmology.h"
#include "Posterior.h"

using namespace std;

std::string cbl::par::DirCosmo;
std::string cbl::par::DirLoc;

// this example model has 4 parameters: A, B, C, D; C and D are free parameters, A is fixed, B is a derived parameter
vector<double> model_function_plus (const vector<double> x, const shared_ptr<void> modelInput, std::vector<double> &parameter)
{
  // the object Cosmology, used in this example to compute Omega_matter
  cbl::cosmology::Cosmology cosm = *static_pointer_cast<cbl::cosmology::Cosmology>(modelInput);

  vector<double> model(x.size(), 0.);
  for (size_t i=0; i<x.size(); ++i)
    model[i] = parameter[2]*x[i]+parameter[3]+parameter[0]*cosm.Omega_matter(); // the model

  parameter[1] = parameter[2]+parameter[3]+parameter[0]; // parameter[1] is a derived parameter

  return model;
}

// this example model has 4 parameters: A, B, C, D; C and D are free parameters, A is fixed, B is a derived parameter
vector<double> model_function_minus (const vector<double> x, const shared_ptr<void> modelInput, std::vector<double> &parameter)
{
  // the object Cosmology, used in this example to compute Omega_matter
  cbl::cosmology::Cosmology cosm = *static_pointer_cast<cbl::cosmology::Cosmology>(modelInput);

  vector<double> model(x.size(), 0.);
  for (size_t i=0; i<x.size(); ++i)
    model[i] = -parameter[2]*x[i]+parameter[3]+parameter[0]*cosm.Omega_matter(); // the model

  parameter[1] = parameter[2]+parameter[3]+parameter[0];  // parameter[1] is a derived parameter

  return model;
}


// Use this function to pass all the model inputs, in this case an object of class cbl::cosmology::Cosmology
shared_ptr<cbl::statistics::Model1D> getModel1D (const cbl::cosmology::Cosmology cosmology, const int model) 
{
  // number of model parameters
  const int nparameters = 4;

  // names of the model parameters
  const vector<string> parNames = {"A", "B", "C", "D"};

  // vector containing the 4 model parameters
  vector<cbl::statistics::ParameterType> parType(nparameters, cbl::statistics::ParameterType::_Base_);
  parType[1] = cbl::statistics::ParameterType::_Derived_;   

  // set the stuff used to construct the model: here an object of class cosmology, just as an example 
  auto ptr_modelInput = make_shared<cbl::cosmology::Cosmology>(cosmology);

  // construct the model
  if (model == 0) {
    const cbl::statistics::Model1D model(&model_function_plus, nparameters, parType, parNames, ptr_modelInput);
    return make_shared<cbl::statistics::Model1D>(model);
  }
  else if (model == 1) {
    const cbl::statistics::Model1D model(&model_function_minus, nparameters, parType, parNames, ptr_modelInput);
    return make_shared<cbl::statistics::Model1D>(model);
  }
  else {
    std::cout << "Illegal option, exiting!" << std::endl;
    exit(0);
  }
}

// Use this function to pass all the model inputs, in this case an object of class cbl::cosmology::Cosmology
shared_ptr<cbl::statistics::Model1D> getModel1D_correlated (const cbl::cosmology::Cosmology cosmology, const int model) 
{
  // number of model parameters
  const int nparameters = 4;

  // names of the model parameters
  const vector<string> parNames = {"A", "B", "C", "D"};

  // vector containing the 4 model parameters
  vector<cbl::statistics::ParameterType> parType(nparameters, cbl::statistics::ParameterType::_Correlated_);  
  parType[0] = cbl::statistics::ParameterType::_Base_;
  parType[1] = cbl::statistics::ParameterType::_Derived_;
  
  // set the stuff used to construct the model: here an object of class cosmology, just as an example 
  auto ptr_modelInput = make_shared<cbl::cosmology::Cosmology>(cosmology);

  // construct the model
  if (model == 0) {
    const cbl::statistics::Model1D model(&model_function_plus, nparameters, parType, parNames, ptr_modelInput);
    return make_shared<cbl::statistics::Model1D>(model);
  }
  else if (model == 1) {
    const cbl::statistics::Model1D model(&model_function_minus, nparameters, parType, parNames, ptr_modelInput);
    return make_shared<cbl::statistics::Model1D>(model);
  }
  else {
    std::cout << "Illegal option, exiting!" << std::endl;
    exit(0);
  }
}
