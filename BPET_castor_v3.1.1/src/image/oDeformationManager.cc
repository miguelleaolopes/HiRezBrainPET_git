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
  \ingroup  image
  \brief    Implementation of class oDeformationManager
*/

#include "oDeformationManager.hh"
#include "sAddonManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn oDeformationManager
  \brief Constructor of oDeformationManager. Simply set all data members to default values.
*/
oDeformationManager::oDeformationManager()
{
  // Image dimensions
  mp_ID = NULL;
  // Options for deformation
  m_options = "";

  // Deformation objects and associated bool
  mp_Deformation = NULL;
  
  m_UseDeformationResp = false;
  m_UseDeformationCard = false;
  m_UseDeformationIPat = false;

  // Variable indicating the number of transformations
  m_nbTransformations = -1;
  
  // Verbosity
  m_verbose = -1;
  // Data mode
  m_dataMode = -1;
  // Not checked yet
  m_checked = false;
  // Not initialized yet
  m_initialized = false;
  
  #ifdef CASTOR_OMP
  // Barrier for multithreaded reconstruction
  mp_deformationRequirement = NULL;
  #endif
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~oDeformationManager
  \brief Destructor of oDeformationManager. Free memory from all allocated tabs.
*/
oDeformationManager::~oDeformationManager() 
{
  // Delete deformation objects
  if (m_initialized)
    if (mp_Deformation!= NULL) delete mp_Deformation;

  #ifdef CASTOR_OMP
  if (mp_deformationRequirement) delete[] mp_deformationRequirement;
  #endif
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      oDeformationManager::SetMotionType
  \param   a_motionType
  \brief   Set the nature of motion correction (Deformation type macro)
*/
void oDeformationManager::SetMotionType(int a_motionType)
{
  switch (a_motionType)
  {
    case DEF_RESP_MOT : m_UseDeformationResp = true; break;
    case DEF_CARD_MOT : m_UseDeformationCard = true; break;
    case DEF_IPAT_MOT : m_UseDeformationIPat = true; break;
    case DEF_DUAL_MOT : {m_UseDeformationResp = true; m_UseDeformationCard = true;} break;
    default : {m_UseDeformationResp = false;
               m_UseDeformationCard = false;
               m_UseDeformationIPat = false;}
  }
}
    
    
    
    
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CheckParameters
  \brief This function is used to check parameters after the latter
         have been all set using Set functions.
  \return 0 if success, positive value otherwise.
*/
int oDeformationManager::CheckParameters()
{
  // Verbose
  if (m_verbose>=VERBOSE_DEBUG_NORMAL) Cout("oDeformationManager::CheckParameters() -> Check parameters before initialization"<< endl);
  // Check image dimensions
  if (mp_ID==NULL)
  {
    Cerr("***** oDeformationManager::CheckParameters() -> No image dimensions provided !" << endl);
    return 1;
  }
  // Check resp gates
  if (m_nbTransformations<0)
  {
    Cerr("***** oDeformationManager::CheckParameters() -> Wrong number of respiratory gates !" << endl);
    return 1;
  }
  // Check verbosity
  if (m_verbose<0)
  {
    Cerr("***** oDeformationManager::CheckParameters() -> Wrong verbosity level provided !" << endl);
    return 1;
  }
  // Check data mode
  if (m_dataMode<0)
  {
    Cerr("***** oDeformationManager::CheckParameters() -> Wrong data mode provided ! (should be =0 (list-mode) or =1 (histogram)" << endl);
    return 1;
  }
  // Normal end
  m_checked = true;
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn Initialize
  \brief Set the flags for the different motion types and instanciate/initialize deformation objects
         through the ParseOptionsAndInitializeDeformations() private function.
  \return 0 if success, positive value otherwise.
*/
int oDeformationManager::Initialize()
{
  // Forbid initialization without check
  if (!m_checked)
  {
    Cerr("***** oDeformationManager::Initialize() -> Must call CheckParameters() before Initialize() !" << endl);
    return 1;
  }

  // Initialize current gate/involuntary motion index
  m_curMotIdx = 0;

  #ifdef CASTOR_OMP
  // initialize flags for managing parallel execution when image deformation is required 
  mp_deformationRequirement = new bool[mp_ID->GetNbThreadsForProjection()];
  for(INTNB th=0 ; th < mp_ID->GetNbThreadsForProjection() ; th++) mp_deformationRequirement[th] = false;
  #endif
  
  // Return if no deformation
  if (!m_UseDeformationCard && !m_UseDeformationIPat && !m_UseDeformationResp)
  {
    m_initialized = true;
    return 0;
  }

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("oDeformationManager::Initialize() -> Initialize deformations" << endl);
  
  // Parse deformation options and initialize them
  if (ParseOptionsAndInitializeDeformations())
  {
    Cerr("***** oDeformationManager::Initialize() -> A problem occurred while parsing deformation options and initializing them !" << endl);
    return 1;
  }

  // Normal end
  m_initialized = true;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ParseOptionsAndInitializeProjectors
  \brief Parse respiratory/cardiac/involuntary patient motion options contained in the previously provided
         strings. This function is called inside the Initialize() function.
  \details Manage the options reading and initialize specific vDeformation
           Options are a string containing first the name of the deformation,
           then either a ':' and a configuration file specific to the deformation
           - or - as many ',' as needed parameters for this deformation.
           Specific pure virtual functions of the vDeformation are used to read parameters and initialize them.
  \todo  Some cleaning if we merge respiratory and cardiac motion objects
  \return 0 if success, positive value otherwise
*/
int oDeformationManager::ParseOptionsAndInitializeDeformations()
{
  if (m_verbose>=VERBOSE_DETAIL) Cout("oDeformationManager::ParseOptionsAndInitializeDeformations ..."<< endl);
    
  string deformation = "";
  string list_options = "";
  string file_options = "";

  // This is for the automatic initialization of the deformations
  typedef vDeformation *(*maker_deformation) ();
  
  // Get deformation's list from addon manager
  std::map <string,maker_deformation> list = sAddonManager::GetInstance()->mp_listOfDeformations;

  size_t colon, comma;
  
  // ---------------------------------------------------------------------------------------------------
  // Manage deformation for respiratory motion
  // ---------------------------------------------------------------------------------------------------
  if (m_UseDeformationResp || m_UseDeformationCard || m_UseDeformationIPat)
  {    
    // ______________________________________________________________________________
    // Get the deformation name in the options and isolate the real deformation's options

    // Search for a colon ":", this indicates that a configuration file is provided after the deformation name
    colon = m_options.find_first_of(":");
    comma = m_options.find_first_of(",");

    // Case 1: we have a colon
    if (colon!=string::npos)
    {
      // Get the deformation name before the colon
      deformation = m_options.substr(0,colon);
      // Get the configuration file after the colon
      file_options = m_options.substr(colon+1);
      // List of options is empty
      list_options = "";
    }
    // Case 2: we have a comma
    else if (comma!=string::npos)
    {
      // Get the deformation name before the first comma
      deformation = m_options.substr(0,comma);
      // Get the list of options after the first comma
      list_options = m_options.substr(comma+1);
      // Configuration file is empty
      file_options = "";
    }
    // Case 3: no colon and no comma (a single deformation name)
    else
    {
      // Get the deformation name
      deformation = m_options;
      // Configuration file is empty
      file_options = "";
      // List of options is empty
      list_options = "";
    }
  
    // Create the deformation
    if (list[deformation]) mp_Deformation = list[deformation]();
    else
    {
      Cerr("***** oDeformationManager::ParseOptionsAndInitializeDeformations() -> Deformation '" << deformation << "' does not exist !" << endl);
      sAddonManager::GetInstance()->ShowHelpDeformation();
      return 1;
    }
    // Set parameters
    mp_Deformation->SetImageDimensionsAndQuantification(mp_ID);
    mp_Deformation->SetVerbose(m_verbose);
    mp_Deformation->SetNbTransformations(m_nbTransformations);
  
    // Provide configuration file if any
    if (file_options!="" && mp_Deformation->ReadAndCheckConfigurationFile(file_options))
    {
      Cerr("***** oDeformationManager::ParseOptionsAndInitializeDeformations() -> A problem occurred while reading and checking respiratory deformation's configuration file !" << endl);
      return 1;
    }
    // Provide options if any
    if (list_options!="" && mp_Deformation->ReadAndCheckOptionsList(list_options))
    {
      Cerr("***** oDeformationManager::ParseOptionsAndInitializeDeformations() -> A problem occurred while parsing and reading respiratory deformation's options !" << endl);
      return 1;
    }
  
    // Check parameters
    if (mp_Deformation->CheckParameters())
    {
      Cerr("***** oDeformationManager::ParseOptionsAndInitializeDeformations() -> A problem occurred while checking respiratory deformation parameters !" << endl);
      return 1;
    }
    // Initialize the deformation
    if (mp_Deformation->Initialize())
    {
      Cerr("***** oDeformationManager::ParseOptionsAndInitializeDeformations() -> A problem occurred while initializing respiratory deformation !" << endl);
      return 1;
    }
    
  }


  // Normal end
  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InstantiateImageForDeformation
  \param oImageSpace* ap_Image : required to call oImageSpace instanciation functions
  \brief If deformation is enabled, ask the Image Space to Instantiate the temporary backward image for deformation
         If reconstruction is in histogram mode, the temporal sensitivity image is instanciated as well
*/
void oDeformationManager::InstantiateImageForDeformation(oImageSpace* ap_Image)
{
  if ( m_UseDeformationResp 
    || m_UseDeformationCard 
    || m_UseDeformationIPat )
  {
    // Verbose
    if(m_verbose>=VERBOSE_DETAIL) Cout("oDeformationManager::InstantiateImageForDeformation() -> Instantiate buffer images"<< endl);
    // Instantiate reference forward/backward images
    ap_Image->InstantiateRefImagesForDeformation();
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn DeallocateImageForDeformation
  \param oImageSpace* ap_Image : required to call oImageSpace deallocation functions
  \brief If deformation is enabled, ask the Image Space to free memory of the temporary backward image for deformation
         If reconstruction is in histogram mode, the temporal sensitivity image is deallocated as well
*/ 
void oDeformationManager::DeallocateImageForDeformation(oImageSpace* ap_Image)
{
  if ( m_UseDeformationResp || m_UseDeformationCard || m_UseDeformationIPat ) 
  {
    // Verbose
    if(m_verbose>=VERBOSE_DETAIL) Cout("oDeformationManager::DeallocateImageForDeformation() -> Deallocate buffer images"<< endl);

    // Free reference forward/backward images
    ap_Image->DeallocateRefImagesForDeformation();
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitImageForDeformation
  \param oImageSpace* ap_Image : required to call oImageSpace initialization functions
  \brief If deformation is enabled, ask the Image Space to initialize the temporary backward image for deformation
         If reconstruction is in histogram mode, the temporal sensitivity image is initialized as well
*/
int oDeformationManager::InitImageForDeformation(oImageSpace* ap_Image)
{
  if ( m_UseDeformationResp || m_UseDeformationCard || m_UseDeformationIPat ) 
  {
    // Reset motion idx
    m_curMotIdx = 0;
    
    if(m_verbose>=VERBOSE_DETAIL) Cout("oDeformationManager::InitImageForDeformation() -> Initialize buffer images"<< endl);

    // Init reference forward/backward images
    ap_Image->InitRefImagesForDeformation();
    
    // For cyclic motion, perform deformation of the forward image
    // on each potential dynamic image
    if (m_UseDeformationResp 
     || m_UseDeformationCard
     || m_UseDeformationIPat)
    {
      if(m_verbose >=VERBOSE_DETAIL) Cout("oDeformationManager::InitImageForDeformation()-> Apply forward image deformation of the forward image(s) for each dynamic image (is any) " << endl);
      
      // Deformation of the forward image for each dynamic image
      for(int fr=0; fr<mp_ID->GetNbTimeBasisFunctions(); fr++)
        for(int rimg=0; rimg<mp_ID->GetNbRespBasisFunctions(); rimg++)
          for(int cimg=0; cimg<mp_ID->GetNbCardBasisFunctions(); cimg++)
            if (mp_Deformation->ApplyDeformations(ap_Image->m4p_refDynForwardImage[ fr ][ rimg ][ cimg ],
                                                  ap_Image->m4p_forwardImage[ fr ][ rimg ][ cimg ], 
                                                  FORWARD_DEFORMATION, 
                                                  // cyclic motion: transformation params #0
                                                  // patient motion-> first transformation params for the frame fr (note: rimg & cimg=0)
                                                  0 + mp_ID->GetPMotionFirstIndexForFrame( fr ) ) 
                                                  )
            {
              Cerr("***** oDeformationManager::InitImageForDeformation()-> An error occurred while performing image deformation during the reconstruction !" << endl);
              return 1;
            }
    }
  }
  
  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ApplyDeformationForSensitivityGeneration
  \param oImageSpace* ap_Image : required to access oImageSpace image matrices
  \param int a_defDirection : direction of the deformation (forward/backward)
  //\param int a_defType : Nature of the motion (Respiratory/Cardiac/Involuntary Patient)
  \param int idx
  \param int fr
  \param int rg
  \param int cg
  \brief Apply deformations during the list-mode sensitivity image generation 
  \details Perform deformation on the forward_image or the backward_image matrices corresponding to the current fr, rg, cg (if any), and depending on the defDirection.
  \todo Some changes required if we merge respiratory/cardiac motion objects
  \todo Check and implement patient motion
  \return 0 if success, positive value otherwise
*/
int oDeformationManager::ApplyDeformationForSensitivityGeneration(oImageSpace* ap_Image, int a_defDirection, int idx, int fr, int rg, int cg) 
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oDeformationManager::ApplyDeformationForSensitivityGeneration() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  if (m_UseDeformationResp)
    idx = rg;
  else if (m_UseDeformationCard)
    idx = cg;
  else if (m_UseDeformationIPat)
    idx += rg; // First motion parameter for the fr + motion image index
  else
    return 0; // no deformation to perform
    

  // --- DEFORMATION MANAGEMENT ---
      
  if (m_verbose>=VERBOSE_DETAIL)
  {
    Cout("oDeformationManager::ApplyDeformationForSensitivityGeneration(): ");
    
    if(a_defDirection==BACKWARD_DEFORMATION)
      Cout("->BACKWARD Deformation of the sensitivity image with transformation parameter "<<idx<<", (respiratory gate "<<rg << " and card gate " <<cg<< ")." <<endl);
    else
      Cout("->FORWARD Deformation of the sensitivity image with transformation parameter "<<idx<<", (respiratory gate "<<rg << " and card gate " <<cg<< ")." <<endl);
  }
  
  if (mp_Deformation->PerformSensitivityDeformation(ap_Image, a_defDirection, idx, fr, rg, cg) )
  {
    Cerr("***** oDeformationManager::ApplyDeformationForSensitivityGeneration()-> An error occurred while performing deformation for respiratory motion correction for the sensitivity image generation !" << endl);
    return 1;
  }
  
  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PerformDeformation
  \param oImageSpace* ap_Image : required to access oImageSpace image matrices
  \brief Apply deformations during reconstruction
  \details Call the eponym function for the deformation object, 
           as well as PerformHistoSensitivityDeformation() if data mode is histogram.
  \todo why the check on frames ?
  \return 0 if success, positive value otherwise
*/
int oDeformationManager::PerformDeformation(oImageSpace* ap_Image)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oDeformationManager::PerformDeformation() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  // Get the deformation index corresponding to the motion 
  int idx = -1;
  
  if (m_UseDeformationResp)
    idx = mp_ID->GetCurrentRespGate(0);
  else if (m_UseDeformationCard)
    idx = mp_ID->GetCurrentCardGate(0);
  else if (m_UseDeformationIPat)
    idx = mp_ID->GetCurrentPMotionIndex(0);
  else
    return 0; // no deformation to perform


  // --- DEFORMATION MANAGEMENT ---

  //int fr = mp_ID->GetCurrentTimeFrame(0); // Done via time_basis_coef
  int rimg = mp_ID->GetCurrentRespImage(0);
  int cimg = mp_ID->GetCurrentCardImage(0);  
  
  // This loop is only required when intrinsic temporal basis functions are used
  // The image matrices are then not independent for each frame, but for each basis functions.
  // In this case, the deformation operations are performed for each image matrix
  // In normal case, the loop is simply ignored as only one frame has time basis coefficient == 1
  // (Technically, it is has using identity intrinsic basis functions)
  for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
  {
    FLTNB time_basis_coef = mp_ID->GetTimeBasisCoefficient(tbf, mp_ID->GetCurrentTimeFrame(0));
    if (time_basis_coef==0.) continue;
    
    // Perform deformation if the motion index has been changed
    if ( m_curMotIdx != (idx+1) )
    {
      if(m_verbose >=VERBOSE_DETAIL) Cout("oDeformationManager::PerformDeformation-> Transformation parameter # " << idx
                                     << ", for resp image " <<rimg
                                     << " and card image " <<cimg<<"." 
                                     << " Frame # " << tbf << endl);
      if(mp_Deformation->PerformDeformation(ap_Image, idx, tbf, rimg, cimg) )
      {
        Cerr("***** oDeformationManager::PerformDeformation()-> An error occurred while performing image deformation during the reconstruction !" << endl);
        return 1;
      }
  
      // Perform deformation of the sensitivity image ((PET) MODE_HISTOGRAM)
      if(m_dataMode == MODE_HISTOGRAM 
      && mp_Deformation->PerformHistoSensitivityDeformation(ap_Image, idx, tbf, rimg, cimg) )
      {
        Cerr("***** oDeformationManager::PerformDeformation()-> An error occurred while performing sensitivity image deformation (histogram mode) during the reconstruction !" << endl);
        return 1;
      }
  
      // Update the current index
      m_curMotIdx = idx;
    }
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ApplyDeformationsToBackwardImage
  \param oImageSpace* ap_Image : required to access oImageSpace image matrices
  \brief Apply final backward deformations on the backward image
  \details Call the eponym function for the deformation object, as well as >ApplyDeformationsToHistoSensitivityImage() if data mode is histogram.
           Then reinitialize the temporary backup deformation images (the backward image, and the sensitivity image if data mode is histogram) 
  \return 0 if success, positive value otherwise
*/
int oDeformationManager::ApplyDeformationsToBackwardImage(oImageSpace* ap_Image)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oDeformationManager::ApplyDeformationsToBackwardImage() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  // Get the deformation index corresponding to the motion 
  int idx = -1;
  
  // Loop on each dynamic frames
  for(int fr=0; fr<mp_ID->GetNbTimeBasisFunctions(); fr++)
  {
    if (m_UseDeformationResp)
      idx = mp_ID->GetCurrentRespGate(0);
    else if (m_UseDeformationCard)
      idx = mp_ID->GetCurrentCardGate(0);
    else if (m_UseDeformationIPat)
      idx = mp_ID->GetPMotionLastIndexForFrame(fr);
    else
      return 0; // no deformation to perform
      
    // --- DEFORMATION MANAGEMENT ---
      
    if(m_verbose >=VERBOSE_DETAIL) Cout("oDeformationManager::ApplyDeformationsToBackwardImage-> Transformation parameter #" << idx<< "." <<endl);
      
    if(mp_Deformation->ApplyDeformationsToBackwardImage(ap_Image, fr, idx) )
    {
      Cerr("***** oDeformationManager::ApplyDeformationsToBackwardImage()-> An error occurred while performing final backward image deformation !" << endl);
      return 1;
    }
      
    // Perform deformation of the sensitivity image ((PET) MODE_HISTOGRAM)
    if(m_dataMode == MODE_HISTOGRAM
    && mp_Deformation->ApplyDeformationsToHistoSensitivityImage(ap_Image, fr, idx) )
    {
      Cerr("***** oDeformationManager::ApplyDeformationsToBackwardImage()-> An error occurred while performing final backward image deformation of the sensitivity image !" << endl);
      return 1;
    }
  }
  
  // Init reference forward/backward images
  ap_Image->InitRefImagesForDeformation(); 

  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn TestDeformationOnImage
      \param   ap_inputImage : input image to deform
      \param   ap_outputImage : image in which the output of the deformation should be recovered
      \param   a_direction : a direction for the deformation to perform (forward or backward)
      \param   a_defIdx : index of the deformation
  \brief Apply deformation specified by arguments on provided input image, for testing purposes
  \return 0 if success, positive value otherwise
*/
int oDeformationManager::TestDeformationOnImage(FLTNB* ap_inputImage, FLTNB* ap_outputImage, int a_direction, int a_defIdx)
{
  if (m_verbose>=VERBOSE_DETAIL) Cout("oDeformationManager::TestDeformationOnImage ..."<< endl);
  
  if(mp_Deformation->ApplyDeformations(ap_inputImage, ap_outputImage, a_direction, a_defIdx) )
  {
    Cerr("***** oDeformationManager::TestDeformationOnImage()-> An error occurred while testing deformation !" << endl);
    return 1;
  }
  
  return 0;
}

#ifdef CASTOR_OMP
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      oDeformationManager::AllThreadsRequireDeformation
  \brief   Do all the threads require an image deformation?
  \return  True if all the threads require a deformation, false otherwise
*/

bool oDeformationManager::AllThreadsRequireDeformation()
{
  for (INTNB th=0 ; th < mp_ID->GetNbThreadsForProjection() ; th++)
    if (!mp_deformationRequirement[th])
      return false;
  return true;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      oDeformationManager::UnsetDeformationRequirements
  \brief   Unset deformation requirements for all the threads
*/

void oDeformationManager::UnsetDeformationRequirements()
{
  for (INTNB th=0 ; th < mp_ID->GetNbThreadsForProjection() ; th++)
    mp_deformationRequirement[th] = false;
}
#endif
