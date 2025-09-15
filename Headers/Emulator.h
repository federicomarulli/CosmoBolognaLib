/********************************************************************
 *  Copyright (C) 2023 by Pierpaolo Nicolosi & Sofia Contarini *
 *  sofia.contarini3@unibo.it * * This program is free software; you
 *  can redistribute it and/or * modify it under the terms of the GNU
 *  General Public License as * published by the Free Software
 *  Foundation; either version 2 of * the License, or (at your option)
 *  any later version.  * * This program is distributed in the hope
 *  that it will be useful, * but WITHOUT ANY WARRANTY; without even
 *  the implied warranty of * MERCHANTABILITY or FITNESS FOR A
 *  PARTICULAR PURPOSE.  See the * GNU General Public License for more
 *  details.  * * You should have received a copy of the GNU General
 *  Public * License along with this program; if not, write to the
 *  Free * Software Foundation, Inc., * 59 Temple Place - Suite 330,
 *  Boston, MA 02111-1307, USA.  *
 ********************************************************************/

/**
 *  @file Headers/Emulator.h
 *
 *  @brief The class Emulator used to handle the parameters necessary
 *  to emulate functions, e.g. computed with CAMB
 *
 *  This file defines the interface of the class Emulator, used to
 *  handle the parameters necessary to emulate functions for cosmology
 *
 *  @authors Pierpaolo Nicolosi, Sofia Contarini
 *
 *  @authors pierpaolo.nicolosi@studio.unibo.it, sofia.contarini3@unibo.it
 */

#ifndef __CAMB_EMU__
#define __CAMB_EMU__

#include "Path.h"
#include "Constants.h"
#include "Func.h"

/**
 *  @brief The global namespace of the <B> \e CosmoBolognaLib </B>
 *  
 *  The \e cbl namespace contains all the main functions and
 *  classes of the CosmoBolognaLib
 */
namespace cbl
{
  /**
   *  @brief The namespace of the functions and classes for the usage 
   *  of emulators
   *  
   *  The \e catalogue namespace contains all the functions and
   *  classes for the usage of emulators
   */
  namespace emulator
  {
    /**
     *  @class Emulator
     *
     *  @brief The class Emulator
     *
     *  This class defines the interface of the class Emulator
     *
     */
    class Emulator
    {

      /// the mean of the taining input parameters
      std::vector<double> m_input_mean;

      /// the standard deviation of the training input parameters
      std::vector<double> m_input_std;

      /// the trained weights
      std::vector<std::vector<std::vector<double>>> m_weights;

      /// the trained bias
      std::vector<std::vector<double>> m_bias;

      /// the trained betas
      std::vector<std::vector<double>> m_betas;

      /// the trained alphas
      std::vector<std::vector<double>> m_alphas;

      /// the mean of the training output model
      std::vector<double> m_output_mean;

      /// the standard deviation of the training output model
      std::vector<double> m_output_std;

      /// the training modes
      std::vector<double> m_modes;

      /**
       *  @brief function to read a matrix from file
       *
       *  @param file_name the name of the files to read
       *
       *  @param n_columns the columns to read
       *
       *  @param path_name the path to stored file
       *
       *  @return a vector of vectors of double
       */
      std::vector<std::vector<double>> m_loadMatrix (const std::string file_name, const int n_columns, const std::string path_name);

      /**
       *  @brief function to initialize the deep 
       *  neural network's layers: input layer, 
       *  hidden layers and output layer.
       *
       *  @return the layers initialized with zeros
       */
      std::vector<std::vector<double>> m_initializeLayers ();

      /**
       *  @brief function normalize the input parameters
       *
       *  @param input_params the input parameters
       *
       *  @return the normalized inputs
       */
      std::vector<double> m_normalizeRepresentation (const std::vector<double> input_params);

      /**
       *  @brief the activation function of each deep neural 
       *  network's node:
       *  
       * \f[ f_{i}(x_{i}) = \left(\beta_{i} + \frac{1 - \beta_{i}}{
       *  1 + e^{-\alpha_{i}x_{i}}}\right)x_{i} \f]
       *
       * where \f$ \alpha_{i} \f$ and \f$ \beta_{i} \f$ are trained
       * deep neuaral network's parameters (see cbl::emulator::
       * Emulator::loadAlphas and cbl::emulator::Emulator::loadBetas),
       * while \f$ x_{i} \f$ is an element of the previous forward
       * passing operation (see
       * cbl::emulator::Emulator::m_layerOperation)
       *
       *  @param activationLayer the layer to activate
       * 
       *  @param alphas \f$ \alpha={\alpha_{1}, \alpha_{2},...} \f$
       *
       *  @param betas \f$ \beta={\beta_{1}, \beta_{2},...} \f$
       *
       *  @return the activated layer
       */
      std::vector<double> m_activactionFunc (const std::vector<double> activationLayer, std::vector<double> alphas, std::vector<double> betas);

      /**
       *  @brief operation function to pass 
       * forward into the deep neural network
       *
       *  @param input_layer the input layer
       *
       *  @param weights the weights
       *
       *  @param bias the bias
       *
       *  @return the activated layer
       */
      std::vector<double> m_layerOperation (const std::vector<double> input_layer, const std::vector<std::vector<double>> weights, const std::vector<double> bias);

      /**
       *  @brief check if the main parameters have been loaded
       *
       */
      void m_checkParameters ();

      /**
       *  @brief function to load the input deep neural network parameters
       *
       *  @param input_mean the mean of the training input 
       *  parameters
       *
       *  @param input_std the standard deviation of the
       *  training input parameters
       *
       *  @param weights the trained weights
       *
       *  @param bias the trained bias
       *
       *  @param alphas the trained alphas
       *
       *  @param betas the trained betas
       *
       *  @param output_mean the mean of the training output 
       *  model
       *
       *  @param output_std the standard deviation of the 
       *  training output model
       *
       *  @param modes the modes of the training output model
       *
       */
      void m_setParameters (const std::vector<double> input_mean, const std::vector<double> input_std, const std::vector<std::vector<std::vector<double>>> weights, const std::vector<std::vector<double>> bias, const std::vector<std::vector<double>> alphas, const std::vector<std::vector<double>> betas, const std::vector<double> output_mean, const std::vector<double> output_std, const std::vector<double> modes);

    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  
       */
      Emulator () = default;

      /**
       *  @brief constructor of emulator with all the deep neural
       *  network (see https://arxiv.org/abs/2106.03846) trained
       *  parameters as input
       *
       *  @param input_mean the mean of the training 
       *  input parameters
       *
       *  @param input_std the standard deviation of 
       *  the training input parameters
       *
       *  @param weights the trained weights
       *
       *  @param bias the trained bias
       *
       *  @param alphas the trained alphas
       *
       *  @param betas the trained betas
       *
       *  @param output_mean the mean of the training 
       *  output model
       *
       *  @param output_std the standard deviation of 
       *  the training output model
       *
       *  @param modes the training modes
       *
       */
      Emulator (const std::vector<double> input_mean, const std::vector<double> input_std, const std::vector<std::vector<std::vector<double>>> weights, const std::vector<std::vector<double>> bias, const std::vector<std::vector<double>> alphas, const std::vector<std::vector<double>> betas, const std::vector<double> output_mean, const std::vector<double> output_std, const std::vector<double> modes);

      /**
       *
       *  @brief function to load all the deep neural network (see
       *  https://arxiv.org/abs/2106.03846) trained parameters at
       *  once, in order: input parameters, weights, bias, alphas,
       *  betas, output model, modes
       *
       *  @param dir_name the name of the directory where the output
       *  files of the training are stored
       * 
       *  @param dir_path the path of the directory where the training
       *  output is stored, by default is
       *  "CosmoBolognaLib/External/Emulator_trained/"
       *
       *  @param file_name the name of the files to read (vector of
       *  strings for each parameter, all with length 1 except for
       *  weights and bias)
       *
       *  @param n_columns the columns to read in the files (vector of
       *  integers for each parameter, all with length 1 except for
       *  weights and bias)
       *
       */
      Emulator (const std::string dir_name, const std::string dir_path = par::defaultString, const std::vector<std::vector<std::string>> file_name = {{"input_specifics"}, {"weights[0]", "weights[1]", "weights[2]", "weights[3]", "weights[4]"}, {"bias", "bias_last"}, {"alphas"}, {"betas"}, {"output_specifics"}, {"modes"}}, const std::vector<std::vector<int>> n_columns = {{4}, {20, 20, 20, 20, 300}, {20, 1}, {20}, {20}, {300}, {1}});

      /**
       *  @brief default destructor
       */
      ~Emulator () = default;

      ///@}
      

      /**
       *  @name Member functions used to get the private members 
       */
      ///@{

      /**
       *  @brief get the training input parameter mean
       *
       *  @return m_input_mean, that is the mean of the training
       *  input parameters
       */
      inline std::vector<double> input_mean () { return m_input_mean; }

      /**
       *  @brief get the training input parameter standard deviation
       *
       *  @return m_input_std, that is the standard deviation
       *  of the training input parameters
       */
      inline std::vector<double> input_std () { return m_input_std; }

      /**
       *  @brief get the trained weights
       *
       *  @return m_weights, that is the trained weights
       */
      inline std::vector<std::vector<std::vector<double>>> weights () { return m_weights; }

      /**
       *  @brief get the trained bias
       *
       *  @return m_bias, that is the trained bias
       */
      inline std::vector<std::vector<double>> bias () { return m_bias; }

      /**
       *  @brief get the trained alphas
       *
       *  @return m_alphas, that is the trained alphas
       */
      inline std::vector<std::vector<double>> alphas () { return m_alphas; }

      /**
       *  @brief get the trained  betas
       *
       *  @return m_betas, that is the trained betas
       */
      inline std::vector<std::vector<double>> betas () { return m_betas; }

      /**
       *  @brief get the training output parameter mean
       *
       *  @return m_output_mean, that is the 
       *  mean of the training output parameters
       */
      inline std::vector<double> output_mean () { return m_output_mean; }

      /**
       *  @brief get the training output parameter standard deviation
       *
       *  @return m_output_std, that is the standard deviation
       *  of the training output parameters
       */
      inline std::vector<double> output_std () { return m_output_std; }

      /**
       *  @brief get the training modes
       *
       *  @return m_modes, the training modes 
       */
      inline std::vector<double> modes () { return m_modes; }

      ///@}

      /**
       *  @name Member functions used to operate on the emulator object
       */
      ///@{

      /**
       *  @brief function to load the input parameter specifics, the
       *  mean and the standard deviation used to train the deep 
       *  neural network
       *
       *  @param file_name the name of the file to read
       *
       *  @param n_columns the columns to read in the file
       *
       *  @param dir_name the name of the directory where the output
       *  files of the training are stored
       *
       *  @param dir_path the path of the directory where the training
       *  output is stored, by default is
       *  "CosmoBolognaLib/External/Emulator_trained/"
       *
       */
      void loadInput (const std::string file_name, const int n_columns, const std::string dir_name, const std::string dir_path = par::defaultString);

      /**
       *  @brief function to load the trained weights of the deep
       *  neural network
       *
       *  @param file_name the name of the files to read (the size of
       *  the vector must be equal to the number of layers of the
       *  neural network)
       *
       *  @param n_columns the columns to read in the files (the size
       *  of the vector must be equal to the number of layers of the
       *  neural network)
       *
       *  @param dir_name the name of the directory where the output
       *  files of the training are stored
       *
       *  @param dir_path the path of the directory where the training
       *  output is stored, by default is
       *  "CosmoBolognaLib/External/Emulator_trained/"
       *
       */
      void loadWeights (const std::vector<std::string> file_name, const std::vector<int> n_columns, const std::string dir_name, const std::string dir_path = par::defaultString);

      /**
       *  @brief function to load the trained bias of the deep neural network
       *
       *  @param file_name the name of the files to read (the size of
       *  the vector must be 2) 
       *
       *  @param n_columns the columns to read in the files (the size of
       *  the vector must be 2)
       *
       *  @param dir_name the name of the directory where the output
       *  files of the training are stored
       *
       *  @param dir_path the path of the directory where the training
       *  output is stored, by default is
       *  "CosmoBolognaLib/External/Emulator_trained/"
       *
       */
      void loadBias (const std::vector<std::string> file_name, const std::vector<int> n_columns, const std::string dir_name, const std::string dir_path = par::defaultString);

      /**
       *  @brief function to load the trained alphas of the deep neural network
       *
       *  @param file_name the name of the file to read
       *
       *  @param n_columns the columns to read in the file
       *
       *  @param dir_name the name of the directory where the output
       *  files of the training are stored
       *
       *  @param dir_path the path of the directory where the training
       *  output is stored by default is
       *  "CosmoBolognaLib/External/Emulator_trained/"
       *
       */
      void loadAlphas (const std::string file_name, const int n_columns, const std::string dir_name, const std::string dir_path = par::defaultString);

      /**
       *  @brief function to load the trained betas of the deep neural network
       *
       *  @param file_name the name of the file to read
       *
       *  @param n_columns the columns to read in the file
       *
       *  @param dir_name the name of the directory where the output
       *  files of the training are stored
       *
       *  @param dir_path the path of the directory where the training
       *  output is stored, by default is
       *  "CosmoBolognaLib/External/Emulator_trained/"
       *
       */
      void loadBetas (std::string file_name, int n_columns, const std::string dir_name, const std::string dir_path = par::defaultString);

      /**
       *  @brief function to load the output model specifics, the mean
       *  and the standard deviation used to train the deep neural 
       *  network
       *
       *  @param file_name the name of the file to read
       *
       *  @param n_columns the columns to read in the file
       *
       *  @param dir_name the name of the directory where the output
       *  files of the training are stored
       *
       *  @param dir_path the path of the directory where the training
       *  output is stored, by default is
       *  "CosmoBolognaLib/External/Emulator_trained/"
       *
       */
      void loadOutput (const std::string file_name, const int n_columns, const std::string dir_name, const std::string dir_path = par::defaultString);

      /**
       *  @brief function to load the modes used to train the deep
       * neural network
       *
       *  @param file_name the name of the file to read
       *
       *  @param n_columns the columns to read in the file
       *
       *  @param dir_name the name of the directory where the output
       *  files of the training are stored
       *
       *  @param dir_path the path of the directory where the training
       *  output is stored, by default is
       *  "CosmoBolognaLib/External/Emulator_trained/"
       *
       */
      void loadModes(const std::string file_name, const int n_columns, const std::string dir_name, const std::string dir_path = par::defaultString);

      /**
       *  @brief function to compute the emulated model 
       *  (2pt correlation function, power spectrum...)
       *
       *  @param cosmologicalParameters the input cosmological
       *  parameters (the same reported in the emulator info file)
       * 
       *  @param output_modes depending on the emulated model,
       *  e.g. the spatial scales for a emualted 2pt correlation
       *  function
       * 
       *  @return the emulated model (2pt correlation function, power
       *  spectrum...)
       *
       *  @warning The loading of the training parameters is required
       *
       */
      std::vector<double> model (const std::vector<double> cosmologicalParameters, const std::vector<double> output_modes);

      /**
       *  @brief function to compute the emulated model 
       *  (2pt correlation function, power spectrum...)
       *
       *  @param cosmologicalParameters the input cosmological
       *  parameters (the same reported in the emulator info file)
       * 
       *  @return a 2D vector containing the scales in the first dimension
       *  and the emulated model (2pt correlation function, power
       *  spectrum...) in the second
       *
       *  @warning The loading of the training parameters is required
       *
       */
      std::vector<std::vector<double>> model (const std::vector<double> cosmologicalParameters);

      ///@}
      
    };
    
  } // namespace emulator
} // namespace cbl

#endif
