/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
 *  federico.marulli3@unibo.it                                      *
 *                                                                  *
 *  This program is free software; you can redistribute it and/or   *
 *  modify it under the terms of the GNU General Public License as  *
 *  published by the Free Software Foundation; either version 2 of  *
 *  the License, or (at your option) any later version.             *
 *                                                                  *
 *  This program is distributed in the hope that it will be useful, *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
 *  GNU General Public License for more details.                    *
 *                                                                  *
 *  You should have received a copy of the GNU General Public       *
 *  License along with this program; if not, write to the Free      *
 *  Software Foundation, Inc.,                                      *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.       *
 ********************************************************************/

/** @file Statistics/Likelihood.cpp
 *
 *  @brief Methods of the class Likelihood 
 *
 *  This file contains the implementation of the methods of the class
 *  Likelihood
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "LikelihoodParameters.h"
#include "Likelihood.h"

using namespace std;

using namespace cbl;

// ============================================================================================


void cbl::statistics::Likelihood::m_set_grid_likelihood_1D (const int npoints, const std::vector<std::vector<double>> parameter_limits, const std::string output_file)
{
  vector<double> pp(npoints), log_ll(npoints);
  double binSize = (parameter_limits[0][1]-parameter_limits[0][0])/(npoints-1);
  for (int i=0; i<npoints; i++){
    pp[i] = parameter_limits[0][0]+binSize*i;
    vector<double> par = {pp[i]};
    log_ll[i] = this->log(par);
  }

  shared_ptr<statistics::STR_likelihood_inputs> likelihood_inputs = static_pointer_cast<statistics::STR_likelihood_inputs>(m_likelihood_inputs);
  likelihood_inputs->interp_function1D = make_shared<glob::FuncGrid>(glob::FuncGrid(pp, log_ll, "Spline"));

  m_log_likelihood_function_grid = LogLikelihood_1D_interpolated;
  m_likelihood_function_grid = [&] (vector<double> &par, const shared_ptr<void> input) { return exp(m_log_likelihood_function_grid(par, input));};

  if(output_file!= par::defaultString){
    ofstream fout(output_file); checkIO(fout, output_file);
    for (int i=0; i<npoints; i++)
      fout << pp[i] << "  " << log_ll[i] << endl;
    
    fout.close();
  }
}


// ============================================================================================


void cbl::statistics::Likelihood::m_set_grid_likelihood_2D (const int npoints, const std::vector<std::vector<double>> parameter_limits, const std::string output_file)
{
  vector<double> pp1(npoints), pp2(npoints);
  vector<vector<double>> log_ll(npoints, vector<double>(npoints, 0));
  double binSize1 = (parameter_limits[0][1]-parameter_limits[0][0])/(npoints-1);
  double binSize2 = (parameter_limits[1][1]-parameter_limits[1][0])/(npoints-1);
  for (int i=0; i<npoints; i++){
    pp1[i] = parameter_limits[0][0]+binSize1*i;
    for (int j=0; j<npoints; j++){
      pp2[j] = parameter_limits[1][0]+binSize2*j;
      vector<double> par = {pp1[i], pp2[j]};
      log_ll[i][j] = this->log(par);
    }
  }

  shared_ptr<statistics::STR_likelihood_inputs> likelihood_inputs = static_pointer_cast<statistics::STR_likelihood_inputs>(m_likelihood_inputs);
  likelihood_inputs->interp_function2D = make_shared<glob::FuncGrid2D>(glob::FuncGrid2D(pp1, pp2, log_ll, "Spline"));

  m_log_likelihood_function_grid = LogLikelihood_2D_interpolated;
  m_likelihood_function_grid = [&] (vector<double> &par, const shared_ptr<void> input) { return exp(m_log_likelihood_function_grid(par, input));};

  if(output_file!= par::defaultString){
    ofstream fout(output_file); checkIO(fout, output_file);
    for (int i=0; i<npoints; i++){
      for (int j=0; j<npoints; j++)
	fout << pp1[i] << "  " << pp2[j] << " " << log_ll[i][j] << endl;
      fout << endl;
    }
    
    fout.close();
  }

}


// ============================================================================================


void cbl::statistics::Likelihood::m_set_grid_likelihood_1D (const std::string input_file)
{
  vector<double> pp, log_ll;

  read_vector(input_file, pp, log_ll);

  shared_ptr<statistics::STR_likelihood_inputs> likelihood_inputs = static_pointer_cast<statistics::STR_likelihood_inputs>(m_likelihood_inputs);
  likelihood_inputs->interp_function1D = make_shared<glob::FuncGrid>(glob::FuncGrid(pp, log_ll, "Spline"));

  m_log_likelihood_function_grid = LogLikelihood_1D_interpolated;
  m_likelihood_function_grid = [&] (vector<double> &par, const shared_ptr<void> input) { return exp(m_log_likelihood_function_grid(par, input));};

}


// ============================================================================================


void cbl::statistics::Likelihood::m_set_grid_likelihood_2D (const std::string input_file)
{
  vector<double> pp1, pp2;
  vector<vector<double>> log_ll;

  read_matrix(input_file, pp1, pp2, log_ll);

   shared_ptr<statistics::STR_likelihood_inputs> likelihood_inputs = static_pointer_cast<statistics::STR_likelihood_inputs>(m_likelihood_inputs);
  likelihood_inputs->interp_function2D = make_shared<glob::FuncGrid2D>(glob::FuncGrid2D(pp1, pp2, log_ll, "Spline"));

  m_log_likelihood_function_grid = LogLikelihood_2D_interpolated;
  m_likelihood_function_grid = [&] (vector<double> &par, const shared_ptr<void> input) { return exp(m_log_likelihood_function_grid(par, input));};
}


// ============================================================================================


cbl::statistics::Likelihood::Likelihood (const std::shared_ptr<data::Data> data, const std::shared_ptr<Model> model, const LikelihoodType likelihood_type, const std::vector<size_t> x_index, const int w_index, const std::shared_ptr<ModelParameters> model_parameters)
{
  set_data(data);
  set_model(model, model_parameters);
  set_function(likelihood_type, x_index, w_index);
}


// ============================================================================================


double cbl::statistics::Likelihood::operator() (std::vector<double> &pp) const
{
  return (m_use_grid) ? m_likelihood_function_grid(pp, m_likelihood_inputs) : m_likelihood_function(pp, m_likelihood_inputs);
}

// ============================================================================================


double cbl::statistics::Likelihood::log (std::vector<double> &pp) const
{
  return (m_use_grid) ? m_log_likelihood_function_grid(pp, m_likelihood_inputs) : m_log_likelihood_function(pp, m_likelihood_inputs);
}


// ============================================================================================


void cbl::statistics::Likelihood::set_data (const shared_ptr<data::Data> data)
{
  m_data = data;
}


// ============================================================================================


void cbl::statistics::Likelihood::set_model (const shared_ptr<Model> model, const shared_ptr<ModelParameters> model_parameters)
{
  switch (model->dimension()) {
    case Dim::_1D_:
      m_model = make_shared<Model1D>(*static_pointer_cast<Model1D>(model));
      break;
    case Dim::_2D_:
      m_model = make_shared<Model2D>(*static_pointer_cast<Model2D>(model));
      break;
    default:
      ErrorCBL("Error in cbl::statistics::Likelihood::set_model() of Likelihood.cpp: dimension shoud be Dim::_1D_ or Dim::_2D_!");
  }

  m_model_parameters = (model_parameters == NULL) ?  make_shared<LikelihoodParameters>(LikelihoodParameters(m_model->parameters()->nparameters(), m_model->parameters()->type(), m_model->parameters()->name())) : model_parameters;
  m_model->set_parameters(m_model_parameters);	
}


// ============================================================================================


void cbl::statistics::Likelihood::unset_grid ()
{
  m_use_grid = false;
}


// ============================================================================================


void cbl::statistics::Likelihood::set_grid (const int npoints, const std::vector<std::vector<double>> parameter_limits, const string file, const bool read)
{
  if (m_likelihood_type==statistics::LikelihoodType::_NotSet_)
    ErrorCBL("Error in cbl::statistics::Likelihood::set_grid() of Likelihood.cpp: the Likelihood function is not set!");

  unsigned int npar = m_model->parameters()->nparameters_free();
  
  if (npar==0)
    ErrorCBL("Error in cbl::statistics::Likelihood::set_grid() of Likelihood.cpp: there is no parameter free to vary");
  if (npar>2)
    ErrorCBL("Error in of Likelihood.cpp: wrong size for the vector of starting parameters");
  if (parameter_limits.size()!=npar)
    ErrorCBL("Error in cbl::statistics::Likelihood::set_grid() of Likelihood.cpp: wrong size for the vector of parameter limits");

  coutCBL << "Computing tabulated likelihood!" << endl;

  if (npar==1) {
    if (read)
      m_set_grid_likelihood_1D (file);
    else
      m_set_grid_likelihood_1D (npoints, parameter_limits, file);
  }
  else if (npar==2) {
    if (read)
      m_set_grid_likelihood_2D (file);
    else
      m_set_grid_likelihood_2D (npoints, parameter_limits, file);
  }

  m_use_grid = true;
  coutCBL << "Done!" << endl;
}

    
// ============================================================================================


void cbl::statistics::Likelihood::set_function (const LikelihoodType likelihood_type, const std::vector<size_t> x_index, const int w_index)
{
  m_x_index = x_index;
  m_w_index = w_index;
  m_likelihood_type = likelihood_type;

  if (m_data->dataType()==data::DataType::_1D_ || m_data->dataType()==data::DataType::_1D_extra_) {

    switch (m_likelihood_type)
    {
      case (LikelihoodType::_Gaussian_Error_):
	
	// check the errors > 0
	for (int i=0; i<m_data->ndata(); i++) {
	  if (m_data->error(i) > 0)
	    continue;
	  else
	    ErrorCBL("Error in cbl::statistics::Likelihood::set_function of Likelihood(). Errors < 0 not allowed!");
	}

	m_log_likelihood_function = &LogLikelihood_Gaussian_1D_error;
	break; 

      case (LikelihoodType::_Gaussian_Covariance_):
	//Also invert covariance!
	m_data->invert_covariance();
	m_log_likelihood_function = &LogLikelihood_Gaussian_1D_covariance;
	break;

      case (LikelihoodType::_Poissonian_):
	m_log_likelihood_function = &LogLikelihood_Poissonian_1D_;
	break; 

      default:
	ErrorCBL("Error in cbl::statistics::Likelihood::set_function() of Likelihood.cpp: type of likelihood not recognized or not yet implemented!");
	break;
    }
  }

  else if (m_data->dataType() == data::DataType::_2D_ || m_data->dataType() == data::DataType::_2D_extra_) {
    switch (m_likelihood_type)
    {
      case (LikelihoodType::_Gaussian_Error_):

	// check the errors > 0
	for (int i=0; i<m_data->xsize(); i++) {
	  for (int j=0; j<m_data->ysize(); j++) {
	    if (m_data->error(i, j) > 0)
	      continue;
	    else
	      ErrorCBL("Error in cbl::statistics::Likelihood::set_function of Likelihood(). Errors < 0 not allowed!");
	  }
	}
	m_log_likelihood_function = &LogLikelihood_Gaussian_2D_error;
	break; 

      case (LikelihoodType::_Poissonian_):
	m_log_likelihood_function = &LogLikelihood_Poissonian_2D_;
	break; 

      default:
	ErrorCBL("Error in cbl::statistics::Likelihood::set_function() of Likelihood.cpp: type of likelihood not recognized or not yet implemented!");
	break;
    }
  }
  else ErrorCBL("Error in cbl::statistics::Likelihood::set_function() of Likelihood.cpp: data type not recognized or not yet implemented!");

  m_likelihood_function = [&] (vector<double> &par, const shared_ptr<void> input) { return exp(m_log_likelihood_function(par, input));};


}


// ============================================================================================


void cbl::statistics::Likelihood::set_function (const Likelihood_function likelihood_function)
{
  m_likelihood_type = LikelihoodType::_UserDefined_;
  m_likelihood_function = likelihood_function;
  m_log_likelihood_function = [&] (vector<double> &par, const shared_ptr<void> input) { return std::log(m_likelihood_function(par, input)); };
}


// ============================================================================================


void cbl::statistics::Likelihood::maximize (const std::vector<double> start, const std::vector<std::vector<double>> parameter_limits, const unsigned int max_iter, const double tol, const double epsilon)
{
  if (m_likelihood_type==statistics::LikelihoodType::_NotSet_)
    ErrorCBL("Error in cbl::statistics::Likelihood::maximize() of Likelihood.cpp: the Likelihood function is not set!");

  vector<double> starting_par;
  vector<vector<double>> limits_par;

  unsigned int npar_free = m_model->parameters()->nparameters_free();
  unsigned int npar = m_model->parameters()->nparameters();

  if (start.size()==npar_free && parameter_limits.size()==npar_free) 
  {
    starting_par=start;
    limits_par = parameter_limits;
  }
  else if (start.size()==npar && parameter_limits.size()==npar)
  {
    for(size_t i=0; i<npar_free; i++){
      starting_par.push_back(start[m_model->parameters()->free_parameters()[i]]);
      limits_par.push_back(parameter_limits[m_model->parameters()->free_parameters()[i]]);
    }
  }
  else
    ErrorCBL("Error in cbl::statistics::Likelihood::maximize() of Likelihood.cpp: Check your inputs!");

  m_likelihood_inputs = make_shared<STR_likelihood_inputs>(STR_likelihood_inputs(m_data, m_model, m_x_index, m_w_index));

  function<double(vector<double> &)> ll = [this](vector<double> & pp) { 
    return -m_log_likelihood_function(pp, m_likelihood_inputs); 
  };

  coutCBL << "Maximizing the likelihood..." << endl;
  vector<double> result = cbl::gsl::GSL_minimize_nD(ll, starting_par, limits_par, max_iter, tol, epsilon);
  coutCBL << "Done!" << endl << endl;

  m_model_parameters->set_bestfit_value(result);
  m_model_parameters->write_bestfit_info();
  coutCBL << "log(Likelihood) = " << this->operator()(result) << endl << endl;
}


// ============================================================================================


void cbl::statistics::Likelihood::write_model (const string output_dir, const string output_file, const std::vector<double> parameters, const std::vector<double> xx, const std::vector<double> yy)
{
  switch (m_model->dimension()) {

    case Dim::_1D_: 
      {
	vector<double> xvec = xx;
	if (xx.size()==0)
	  xvec = m_data->xx();

	m_model->write(output_dir, output_file, xvec, parameters);
      }

      break;
    case Dim::_2D_: 
      {
	vector<double> xvec = xx, yvec = yy;
	if (xx.size()==0)
	  xvec = m_data->xx();
	if (yy.size()==0)
	  yvec = m_data->yy();

	m_model->write(output_dir, output_file, xvec, yvec, parameters);
      }

      break;
    default:
      ErrorCBL("Error in cbl::statistics::Likelihood::write_model() of set_model.cpp: dimension shoud be Dim::_1D_ or Dim::_2D_!");
  }
}
