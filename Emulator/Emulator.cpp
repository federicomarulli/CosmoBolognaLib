/********************************************************************
 *  Copyright (C) 2023 by Pierpaolo Nicolosi & Sofia Contarini      *
 *  sofia.contarini3@unibo.it                                       *
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

/**
 *  @file Emulator/Emulator.cpp
 *
 *  @brief The class Emulator used to handle the parameters necessary
 *  to emulate functions comupted with CAMB
 *
 *  This file defines the interface of the class Emulator, used to handle
 *  the parameters necessary to emulate functions comupted with CAMB
 *
 *  @authors Pierpaolo Nicolosi, Sofia Contarini
 *
 *  @authors pierpaolo.nicolosi@studio.unibo.it, sofia.contarini3@unibo.it
 */

#include "Emulator.h"

using namespace std;
using namespace cbl;
using namespace emulator;


// ============================================================================


cbl::emulator::Emulator::Emulator (const std::vector<double> input_mean, const std::vector<double> input_std, const std::vector<std::vector<std::vector<double>>> weights, const std::vector<std::vector<double>> bias, const std::vector<std::vector<double>> alphas, const std::vector<std::vector<double>> betas, const std::vector<double> output_mean, const std::vector<double> output_std, const std::vector<double> modes)
{
  m_setParameters(input_mean, input_std, weights, bias, alphas, betas, output_mean, output_std, modes);
}


// ============================================================================


cbl::emulator::Emulator::Emulator (const std::string dir_name, const std::string dir_path, const std::vector<std::vector<std::string>> file_name, const std::vector<std::vector<int>> n_columns)
{
  if (file_name.size() == 7 && n_columns.size() == 7)
    {
      string dirInput;

      if (dir_path == "NULL")
	{
	  Path path;
	  dirInput = path.DirCosmo() + "/External/Emulator_trained/";
	}

      else
	dirInput = dir_path;

      char ch;
      string input = dirInput + dir_name + "/info";
      ifstream fin(input.c_str());
      checkIO(fin, input);
      cout << par::col_green << endl;
      while (fin >> noskipws >> ch)
	cout << ch;
      cout << par::col_default << endl;

      loadInput(file_name[0][0], n_columns[0][0], dir_name, dir_path);
      loadWeights(file_name[1], n_columns[1], dir_name, dir_path);
      loadBias(file_name[2], n_columns[2], dir_name, dir_path);
      loadAlphas(file_name[3][0], n_columns[3][0], dir_name, dir_path);
      loadBetas(file_name[4][0], n_columns[4][0], dir_name, dir_path);
      loadOutput(file_name[5][0], n_columns[5][0], dir_name, dir_path);
      loadModes(file_name[6][0], n_columns[6][0], dir_name, dir_path);
    }

  else
    ErrorCBL("The first size of file_name and n_columns must be 7!", "Emulator", "Emulator.h");
}


// ============================================================================


void cbl::emulator::Emulator::m_setParameters (const std::vector<double> input_mean, const std::vector<double> input_std, const std::vector<std::vector<std::vector<double>>> weights, const std::vector<std::vector<double>> bias, const std::vector<std::vector<double>> alphas, const std::vector<std::vector<double>> betas, const std::vector<double> output_mean, const std::vector<double> output_std, const std::vector<double> modes)
{
  m_input_mean = input_mean;
  m_input_std = input_std;
  m_weights = weights;
  m_bias = bias;
  m_alphas = alphas;
  m_betas = betas;
  m_output_mean = output_mean;
  m_output_std = output_std;
  m_modes = modes;
}


// ============================================================================


void cbl::emulator::Emulator::m_checkParameters ()
{
  if (!m_input_mean.size())
    ErrorCBL("Load the input_mean parameter!", "m_checkParameters", "Emulator.h");

  else if (!m_input_std.size())
    ErrorCBL("Load the input_std parameter!", "m_checkParameters", "Emulator.h");

  else if (!m_weights.size())
    ErrorCBL("Load the weight parameter!", "m_checkParameters", "Emulator.h");

  else if (!m_bias.size())
    ErrorCBL("Load the bias parameter!", "m_checkParameters", "Emulator.h");

  else if (!m_betas.size())
    ErrorCBL("Load the betas parameter!", "m_checkParameters", "Emulator.h");

  else if (!m_alphas.size())
    ErrorCBL("Load the alphas parameter!", "m_checkParameters", "Emulator.h");

  else if (!m_output_mean.size())
    ErrorCBL("Load the output_mean parameter!", "m_checkParameters", "Emulator.h");

  else if (!m_output_std.size())
    ErrorCBL("Load the output_std parameter!", "m_checkParameters", "Emulator.h");

  else if (!m_modes.size())
    ErrorCBL("Load the modes parameter!", "m_checkParameters", "Emulator.h");
}


// ============================================================================


std::vector<std::vector<double>> cbl::emulator::Emulator::m_loadMatrix (const std::string file_name, const int n_columns, const std::string path_name)
{
  vector<int> column_data = linear_bin_vector(n_columns, 1, n_columns);
  vector<vector<double>> matrix = read_file(file_name, path_name, column_data);
  vector<vector<double>> output_matrix = transpose(matrix);

  return output_matrix;
}


// ============================================================================


std::vector<std::vector<double>> cbl::emulator::Emulator::m_initializeLayers ()
{
  vector<vector<double>> output(m_weights.size());

  for (size_t i = 0; i < output.size(); i++)
    output[i] = vector<double>(m_weights[i][0].size(), 0.);

  return output;
}


// ============================================================================


std::vector<double> cbl::emulator::Emulator::m_normalizeRepresentation (const std::vector<double> input_params)
{
  vector<double> output(input_params.size(), 0.);

  for (size_t j = 0; j < input_params.size(); j++)
    output[j] = (input_params[j] - m_input_mean[j]) / m_input_std[j];

  return output;
}


// ============================================================================


std::vector<double> cbl::emulator::Emulator::m_activactionFunc (const std::vector<double> activationLayer, std::vector<double> alphas, std::vector<double> betas)
{
  vector<double> output(activationLayer.size(), 0.);

  for (size_t j = 0; j < output.size(); j++)
    output[j] = (betas[j] + (1. - betas[j]) / (1. + exp(-alphas[j] * activationLayer[j]))) * activationLayer[j];

  return output;
}


// ============================================================================


std::vector<double> cbl::emulator::Emulator::m_layerOperation (const std::vector<double> input_layer, const std::vector<std::vector<double>> weights, const std::vector<double> bias)
{
  vector<double> output_layer(weights[0].size(), 0.);

  for (size_t k = 0; k < weights.size(); k++)
    for (size_t j = 0; j < output_layer.size(); j++)
      output_layer[j] += input_layer[k] * weights[k][j];

  for (size_t i = 0; i < bias.size(); i++)
    output_layer[i] += bias[i];

  return output_layer;
}


// ============================================================================


void cbl::emulator::Emulator::loadInput (const std::string file_name, const int n_columns, const std::string dir_name, const std::string dir_path)
{
  string dirInput;

  if (dir_path == "NULL")
    {
      Path path;
      dirInput = path.DirCosmo() + "/External/Emulator_trained/";
    }

  else
    dirInput = dir_path;

  vector<vector<double>> input_parameters = m_loadMatrix(file_name, n_columns, dirInput + dir_name + "/");

  if (input_parameters.size() != 2)
    ErrorCBL("Check your inputs!", "loadInput", "Emulator.h");

  m_input_mean = input_parameters[0];
  m_input_std = input_parameters[1];
}


// ============================================================================


void cbl::emulator::Emulator::loadWeights (const std::vector<std::string> file_name, const std::vector<int> n_columns, const std::string dir_name, const std::string dir_path)
{
  string dirInput;

  if (dir_path == "NULL")
    {
      Path path;
      dirInput = path.DirCosmo() + "/External/Emulator_trained/";
    }

  else
    dirInput = dir_path;

  if (file_name.size() != n_columns.size())
    ErrorCBL("The file_name and n_columns must have the same size!", "loadWeights", "Emulator.h");

  size_t layers_number = file_name.size();
  vector<vector<vector<double>>> weights(layers_number);

  for (size_t i = 0; i < layers_number; i++)
    weights[i] = m_loadMatrix(file_name[i], n_columns[i], dirInput + dir_name + "/");

  m_weights = weights;
}


// ============================================================================


void cbl::emulator::Emulator::loadBias (const std::vector<std::string> file_name, const std::vector<int> n_columns, const std::string dir_name, const std::string dir_path)
{
  string dirInput;

  if (dir_path == "NULL")
    {
      Path path;
      dirInput = path.DirCosmo() + "/External/Emulator_trained/";
    }

  else
    dirInput = dir_path;

  if (m_weights.size() == 0)
    ErrorCBL("You must load the weights before the bias parameter!", "loadWeights", "Emulator.h");

  if (file_name.size() == 2 && n_columns.size() == 2)
    {

      vector<vector<double>> bias(m_weights.size());

      vector<vector<double>> bias_firsts = m_loadMatrix(file_name[0], n_columns[0], dirInput + dir_name + "/");

      for (size_t i = 0; i < bias_firsts.size(); i++)
	bias[i] = bias_firsts[i];

      bias[bias.size() - 1] = read_file(file_name[1], dirInput + dir_name + "/", {n_columns[1]})[0];

      m_bias = bias;
    }

  else
    ErrorCBL("The file_name and n_columns must have both size 2", "loadBias", "Emulator.h");
}


// ============================================================================


void cbl::emulator::Emulator::loadAlphas (const std::string file_name, const int n_columns, const std::string dir_name, const std::string dir_path)
{
  string dirInput;

  if (dir_path == "NULL")
    {
      Path path;
      dirInput = path.DirCosmo() + "/External/Emulator_trained/";
    }

  else
    dirInput = dir_path;

  vector<vector<double>> alphas = m_loadMatrix(file_name, n_columns, dirInput + dir_name + "/");

  m_alphas = alphas;
}


// ============================================================================


void cbl::emulator::Emulator::loadBetas (const std::string file_name, const int n_columns, const std::string dir_name, const std::string dir_path)
{
  string dirInput;

  if (dir_path == "NULL")
    {
      Path path;
      dirInput = path.DirCosmo() + "/External/Emulator_trained/";
    }

  else
    dirInput = dir_path;

  vector<vector<double>> betas = m_loadMatrix(file_name, n_columns, dirInput + dir_name + "/");

  m_betas = betas;
}


// ============================================================================


void cbl::emulator::Emulator::loadOutput (const std::string file_name, const int n_columns, const std::string dir_name, const std::string dir_path)
{
  string dirInput;

  if (dir_path == "NULL")
    {
      Path path;
      dirInput = path.DirCosmo() + "/External/Emulator_trained/";
    }

  else
    dirInput = dir_path;

  vector<vector<double>> output_parameters = m_loadMatrix(file_name, n_columns, dirInput + dir_name + "/");

  if (output_parameters.size() != 2)
    ErrorCBL("Check your inputs!", "loadOutput", "Emulator.h");

  m_output_mean = output_parameters[0];
  m_output_std = output_parameters[1];
}


// ============================================================================


void cbl::emulator::Emulator::loadModes (const std::string file_name, const int n_columns, const std::string dir_name, const std::string dir_path)
{
  string dirInput;

  if (dir_path == "NULL")
    {
      Path path;
      dirInput = path.DirCosmo() + "/External/Emulator_trained/";
    }

  else
    dirInput = dir_path;

  vector<double> modes = read_file(file_name, dirInput + dir_name + "/", {n_columns})[0];

  m_modes = modes;
}


// ============================================================================


std::vector<double> cbl::emulator::Emulator::model (const std::vector<double> cosmologicalParameters, const std::vector<double> output_modes)
{
  m_checkParameters();

  std::vector<double> input_params = cosmologicalParameters;

  vector<double> Z_norm_input = m_normalizeRepresentation(input_params);

  vector<vector<double>> layers = m_initializeLayers();
  vector<vector<double>> act(m_betas.size(), vector<double>(m_betas[0].size(), 0.));

  act[0] = m_layerOperation(Z_norm_input, m_weights[0], m_bias[0]);
  layers[0] = m_activactionFunc(act[0], m_alphas[0], m_betas[0]);

  for (size_t i = 1; i < act.size(); i++)
    {
      act[i] = m_layerOperation(layers[i - 1], m_weights[i], m_bias[i]);
      layers[i] = m_activactionFunc(act[i], m_alphas[i], m_betas[i]);
    }

  layers[layers.size() - 1] = m_layerOperation(layers[layers.size() - 2], m_weights[layers.size() - 1], m_bias[layers.size() - 1]);
  for (size_t j = 0; j < m_output_std.size(); j++)
    layers[layers.size() - 1][j] = layers[layers.size() - 1][j] * m_output_std[j] + m_output_mean[j];

  std::vector<double> func = layers[layers.size() - 1];
  std::vector<double> interp_func(output_modes.size(), 0.);

  for (size_t i = 0; i < output_modes.size(); i++)
    interp_func[i] = (interpolated(output_modes[i], m_modes, func, "Linear"));

  return interp_func;
}

// ============================================================================

std::vector<std::vector<double>> cbl::emulator::Emulator::model (const std::vector<double> cosmologicalParameters)
{
  m_checkParameters();

  std::vector<double> input_params = cosmologicalParameters;

  vector<double> Z_norm_input = m_normalizeRepresentation(input_params);

  vector<vector<double>> layers = m_initializeLayers();
  vector<vector<double>> act(m_betas.size(), vector<double>(m_betas[0].size(), 0.));

  act[0] = m_layerOperation(Z_norm_input, m_weights[0], m_bias[0]);
  layers[0] = m_activactionFunc(act[0], m_alphas[0], m_betas[0]);

  for (size_t i=1; i<act.size(); i++) {
    act[i] = m_layerOperation(layers[i - 1], m_weights[i], m_bias[i]);
    layers[i] = m_activactionFunc(act[i], m_alphas[i], m_betas[i]);
  }

  layers[layers.size() - 1] = m_layerOperation(layers[layers.size() - 2], m_weights[layers.size() - 1], m_bias[layers.size() - 1]);
  for (size_t j=0; j<m_output_std.size(); j++)
    layers[layers.size() - 1][j] = layers[layers.size() - 1][j] * m_output_std[j] + m_output_mean[j];

  std::vector<double> func = layers[layers.size() - 1];
  std::vector<std::vector<double>> output_comb(2, std::vector<double>(m_modes.size()));

  output_comb[0] = m_modes;
  output_comb[1] = func;

  return output_comb;
}
