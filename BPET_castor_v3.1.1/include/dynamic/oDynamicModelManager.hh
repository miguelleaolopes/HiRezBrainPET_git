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
  \ingroup  dynamic
  \brief    Declaration of class oDynamicModelManager
*/

#ifndef ODYNAMICMODELMANAGER_HH
#define ODYNAMICMODELMANAGER_HH 1

#include "gVariables.hh"
#include "vDynamicModel.hh"
#include "vDataFile.hh"

/*!
  \class   oDynamicModelManager
  \brief   This class is designed to manage the use of dynamic model in the reconstruction
  \details As each manager class, it is created in the main program, all parameters are
           then set, checked, and the manager is initialized. \n
           The manager is then used by the algorithm itself, where the function
           ApplyDynamicModel() is called at the end of the subset loop
*/
class oDynamicModelManager
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      oDynamicModelManager::oDynamicModelManager
      \brief   Constructor of oDynamicModelManager. Simply set all data members to default values.
    */
    oDynamicModelManager();
    /*!
      \fn      oDynamicModelManager::~oDynamicModelManager
      \brief   Destructor of oDynamicModelManager. Free memory from all allocated tabs.
    */
    ~oDynamicModelManager();

  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      oDynamicModelManager::CheckParameters
      \brief   This function is used to check parameters after the latter
               have been all set using Set functions.
      \return  0 if success, positive value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      oDynamicModelManager::Initialize
      \brief   Set the dynamic model flag and instanciate/initialize model objects
               through the ParseOptionsAndInitializeDeformations() private function.
      \return  0 if success, positive value otherwise.
    */
    int Initialize();
    /*!
      \fn      oDynamicModelManager::oDynamicModelManager::ApplyDynamicModel
      \param   ap_ImageS : pointer to the ImageSpace
      \param   a_ite : index of the actual iteration
      \param   a_sset : index of the actual subset
      \brief   Call EstimateModelParameters() and FitModel() functions of 
               the dynamic model object is 'm_useModel' is on.
      \return  0 if success, positive value otherwise
    */
    int ApplyDynamicModel(oImageSpace* ap_ImageS, int a_iteration, int a_subset);
    /*!
      \fn      oDynamicModelManager::SaveParametricImages
      \param   a_iteration : current iteration index
      \param   a_subset : current number of subsets (or -1 by default)
      \brief   Call SaveParametricImages() function of the dynamic model object is 
              'm_useModel' is on, in order to save any parameter image
      \return  0 if success, positive value otherwise
    */
    int SaveParametricImages(int a_iteration, int a_subset = -1);
    /*!
      \fn      oDynamicModelManager::SetDiagonalBasisFunctions
      \brief   Use SetDiagonalBasisFunctions() when user has not provided a dynamic-model option
               to populate the basis functions in ImageDimensionsAndQuantification for simple.
      \return  0 if success, positive value otherwise
    */


  // -----------------------------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      oDynamicModelManager::SetVerbose
      \param   a_verboseLevel
      \brief   Set the verbose level
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;}
    /*!
      \fn      oDynamicModelManager::SetImageDimensionsAndQuantification
      \param   ap_ImageDimensionsAndQuantification
      \brief   Set the image dimensions in use
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification) 
                                                   {mp_ID = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      oDynamicModelManager::SetOptions
      \param   a_options
      \brief   Set the respiratory motion options contained in the provided string
    */
    inline void SetOptions(const string& a_options)
           {m_options = a_options;}
    /*!
      \fn      oDynamicModelManager::GetUseModel
      \brief   Indicate if the use of a dynamic model is enabled
      \return  true if enabled, false otherwise
    */       
    inline bool GetUseModel()
           {return m_useModel;}
    /*!
      \fn      oDynamicModelManager::SetUseModelInReconstruction
      \brief   Set flag to indicate if the dynamic model
               is used in the tomographic reconstruction
      \return  true if enabled, false otherwise
    */
    inline void SetUseModelInReconstruction(bool a_useModelInReconstruction)
           {m_useModelInReconstruction=a_useModelInReconstruction;}
    /*!
      \fn      oDynamicModelManager::SetDiagonalTimeBasisFunctions
      \brief   Set diagonal Time Basis Functions for regular frame-by-frame reconstruction
      \return  true if enabled, false otherwise
    */
    int SetDiagonalTimeBasisFunctions();
    /*!
      \fn      oDynamicModelManager::SetDiagonalRespBasisFunctions
      \brief   Set diagonal Respiratory Basis Functions for regular gated reconstruction
      \return  true if enabled, false otherwise
    */
    int SetDiagonalRespBasisFunctions();
    /*!
      \fn      oDynamicModelManager::SetDiagonalCardBasisFunctions
      \brief   Set diagonal Cardiac Basis Functions for regular gated reconstruction
      \return  true if enabled, false otherwise
    */
    int SetDiagonalCardBasisFunctions();

  // -----------------------------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      ParseOptionsAndInitializeModel
      \brief   Parse dynamic model options contained in the previously provided
               strings. This function is called inside the Initialize() function.
      \details Manage the options reading and initialize specific vDynamicModel \n
               Options are a string containing first the name of the model,
               then either a ':' and a configuration file specific to the model
               - or - as many ',' as needed parameters for this model. \n
               Specific pure virtual functions of the vDynamicModel are used to read parameters and initialize them.
      \return  0 if success, positive value otherwise
    */
    int ParseOptionsAndInitializeModel();


  // -----------------------------------------------------------------------------------------
  // Data members
  private:

    oImageDimensionsAndQuantification* mp_ID; /*!< Pointer to the oImageDimensionsAndQuantification object in use */
    string m_options;                         /*!< The string containing options for the dynamic model */
    vDynamicModel* mp_DynamicModel;           /*!< Dynamic model object */
    // Variables related to basis functions for the ImageDimensionsAndQuantification object
    int m_nbTimeBasisFunctions;               /*!< Number of Time Basis Functions */
    int m_nbRespBasisFunctions;               /*!< Number of Respiratory Basis Functions */
    int m_nbCardBasisFunctions;               /*!< Number of Cardiac Basis Functions */
    FLTNB** m2p_timeBasisFunctions;           /*!< The table of time basis functions coefficients */
    FLTNB** m2p_respBasisFunctions;           /*!< The table of respiratory basis functions coefficients */
    FLTNB** m2p_cardBasisFunctions;           /*!< The table of cardiac basis functions coefficients */


    bool m_useModel;                          /*!< Flag indicating that use of a dynamic model is enabled or not*/
    bool m_useModelInReconstruction;          /*!< Flag indicating if the model is used with castor-recon or with imageDynamicTools */
    int m_verbose;                            /*!< The verbose level */
    bool m_checked;                           /*!< Boolean indicating whether the parameters were checked or not */
    bool m_initialized;                       /*!< Boolean indicating whether the manager was initialized or not */
};


#endif
