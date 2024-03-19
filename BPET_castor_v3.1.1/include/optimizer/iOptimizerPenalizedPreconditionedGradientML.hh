/*
This file is part of CASToR.

    CASToR is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    CASToR is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along with
    CASToR (in file GNU_GPL.TXT). If not, see <http://www.gnu.org/licenses/>.

Copyright 2017-2021 all CASToR contributors listed below:

    --> Didier BENOIT, Claude COMTAT, Marina FILIPOVIC, Thibaut MERLIN, Mael MILLARDET, Simon STUTE, Valentin VIELZEUF, Zacharias CHALAMPALAKIS

This is CASToR version 3.1.1.
*/

/*!
  \file
  \ingroup  optimizer
  \brief    Declaration of class iOptimizerPenalizedPreconditionedGradient
*/

#ifndef IOPTIMIZERPENALIZEDPRECONDITIONEDGRADIENTML_HH
#define IOPTIMIZERPENALIZEDPRECONDITIONEDGRADIENTML_HH 1

#include "gVariables.hh"
#include "sAddonManager.hh"
#include "vOptimizer.hh"


/*!
  \class   iOptimizerPenalizedPreconditionedGradientML
  \brief   This class implements the Preconditioned Gradient MAP algorithm
  \details This class inherits from vOptimizer and implements the PreconditionedGradientMAP algorithm from Nuyts et al, 2002
*/
class iOptimizerPenalizedPreconditionedGradientML : public vOptimizer
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public iOptimizerPenalizedPreconditionedGradientML::iOptimizerPenalizedPreconditionedGradientML()
      \brief   The constructor of iOptimizerPenalizedPreconditionedGradientML
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    iOptimizerPenalizedPreconditionedGradientML();
    /*!
      \fn      public iOptimizerPenalizedPreconditionedGradientML::~iOptimizerPenalizedPreconditionedGradientML()
      \brief   The destructor of iOptimizerPenalizedPreconditionedGradientML
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
    */
    ~iOptimizerPenalizedPreconditionedGradientML();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameter and do not add semi-column at the end of the line)
    FUNCTION_OPTIMIZER(iOptimizerPenalizedPreconditionedGradientML)
    /*!
      \fn      public int iOptimizerPenalizedPreconditionedGradientML::ReadConfigurationFile()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to the child optimizer, from
               a configuration file. It is the implementation of the pure virtual function inherited
               from the abstract class vOptimizer. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadConfigurationFile(const string& a_configurationFile);
    /*!
      \fn      public int iOptimizerPenalizedPreconditionedGradientML::ReadOptionsList()
      \param   const string& a_optionsList
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to the child optimizer, from
               a list of options. It is the implementation of the pure virtual function inherited
               from the abstract class vOptimizer. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadOptionsList(const string& a_optionsList);

  // -------------------------------------------------------------------
  // Private member functions (virtual in vOptimizer)
  private:
  
    /*!
      \fn      private int iOptimizerPenalizedPreconditionedGradientML::PreImageUpdateSpecificStep()
      \brief   A private function used to compute the penalty term of the PreconditionedGradientMAP algorithm
      \details This function implements the virtual eponym function of vOptimizer.
               It computes the penalty terms of the algorithm.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int PreImageUpdateSpecificStep();

  // -------------------------------------------------------------------
  // Private member functions (pure virtual in vOptimizer)
  private:
    /*!
      \fn      private void iOptimizerPenalizedPreconditionedGradientML::ShowHelpSpecific()
      \brief   A function used to show help about the child optimizer
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the optimizer, and how to set them through the use
               of a configuration file or a list of options. It is the implementation of the pure
               virtual function inherited from the abstract class vOptimizer. It is called by the
               public ShowHelp() function.
    */
    void ShowHelpSpecific();
    /*!
      \fn      private int iOptimizerPenalizedPreconditionedGradientML::CheckSpecificParameters()
      \brief   A private function used to check the parameters settings specific to the child optimizer
      \details This function is used to check that all parameters specific to the optimizer are correctly set
               within allowed values. It is called by the CheckParameters() function of the mother class.
               It is the implementation of the pure virtual function inherited from the abstract mother
               class vOptimizer.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      private int iOptimizerPenalizedPreconditionedGradientML::InitializeSpecific()
      \brief   This function is used to initialize specific stuff to the child optimizer.
      \details It is called by the public Initialize() function from the mother.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int InitializeSpecific();
    /*!
      \fn      private int iOptimizerPenalizedPreconditionedGradientML::SensitivitySpecificOperations()
      \param   FLTNB a_data
      \param   FLTNB a_forwardModel
      \param   FLTNB* ap_weight
      \param   FLTNB a_multiplicativeCorrections
      \param   FLTNB a_additiveCorrections
      \param   FLTNB a_quantificationFactor
      \param   oProjectionLine* ap_Line
      \brief   This function compute the weight associated to the provided event (for sensitivity computation)
      \details It is the implementation of the pure virtual function from vOptimizer. The result is
               put at ap_weight location.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                       FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                       FLTNB a_quantificationFactor, oProjectionLine* ap_Line );
    /*!
      \fn      private int iOptimizerPenalizedPreconditionedGradientML::DataSpaceSpecificOperations()
      \param   FLTNB a_data
      \param   FLTNB a_forwardModel
      \param   FLTNB* ap_backwardValues
      \param   FLTNB a_multiplicativeCorrections
      \param   FLTNB a_additiveCorrections
      \param   FLTNB a_quantificationFactor
      \param   oProjectionLine* ap_Line
      \brief   This function performs the data space operations specific to the optimizer (computes the values
               to be backprojected)
      \details It is the implementation of the pure virtual function from vOptimizer. The results to be
               backprojected is put at ap_backwardValues location.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                     FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                     FLTNB a_quantificationFactor, oProjectionLine* ap_Line );
    /*!
      \fn      private int iOptimizerPenalizedPreconditionedGradientML::ImageSpaceSpecificOperations()
      \param   FLTNB a_currentImageValue
      \param   FLTNB* ap_newImageValue
      \param   FLTNB a_sensitivity
      \param   FLTNB* ap_correctionValues
      \param   INTNB a_voxel
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \brief   This function perform the image update step specific to the optimizer
      \details It is the implementation of the pure virtual function from vOptimizer. The new image value is
               put at the ap_newImageValue location.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                      FLTNB a_sensitivity, FLTNB* ap_correctionValues,
                                      INTNB a_voxel, int a_tbf = -1, int a_rbf = -1, int a_cbf = -1 );


  // -------------------------------------------------------------------
  // Data members
  private:
    FLTNB m_dataSpaceDenominatorThreshold;      /*!< Threshold applied to the denominator of the data space operation to avoid 0-divisions or too high ratios */
    FLTNB m_minimumImageUpdateFactor;           /*!< Minimum allowed image update factor, useful to avoid voxels trapped in 0 value (null or negative values mean no restriction) */
    FLTNB m_maximumImageUpdateFactor;           /*!< Maximum allowed image update factor (null or negative values mean no restriction) */
    FLTNB**** m4p_firstDerivativePenaltyImage;  /*!< Image containing the first order derivative penalty terms of the algorithm */
    FLTNB**** m4p_secondDerivativePenaltyImage; /*!< Image containing the second order derivative penalty terms of the algorithm */
};


// Class for automatic insertion (set here the visible optimizer's name as the first parameter,
// put the class name as the second parameter and do NOT add semi-colon at the end of the line)
CLASS_OPTIMIZER(PPGML,iOptimizerPenalizedPreconditionedGradientML)

#endif

