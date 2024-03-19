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
  \brief    Implementation of class oImageDimensionsAndQuantification
*/

#include "oImageDimensionsAndQuantification.hh"
#include "vDataFile.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oImageDimensionsAndQuantification::oImageDimensionsAndQuantification()
{
  // Set all members to default values
  m_nbVoxX = -1;
  m_nbVoxY = -1;
  m_nbVoxZ = -1;
  m_nbVoxXY = -1;
  m_nbVoxXYZ = -1;
  m_voxSizeX = -1.;
  m_voxSizeY = -1.;
  m_voxSizeZ = -1.;
  m_fovSizeX = -1.;
  m_fovSizeY = -1.;
  m_fovSizeZ = -1.;
  m_fovOutPercent = 0.;
  m_flipOutX = false;
  m_flipOutY = false;
  m_flipOutZ = false;
  m_offsetX = 0.;
  m_offsetY = 0.;
  m_offsetZ = 0.;
  m_frameList = "";
  m_nbFramesToSkip = 0;
  m_nbTimeFrames = 1;
  m_nbTimeBasisFunctions = -1;
  m2p_timeBasisFunctions = NULL;
  m2p_frameDurationsInMs = NULL;
  m2p_frameTimeStartInMs = NULL;
  m2p_frameTimeStopInMs = NULL;
  m_timeStaticFlag = false;
  m3p_quantificationFactors = NULL;
  m_ignoredCorrectionsList = "";
  m_ignoreAttnCorrectionFlag = false;
  m_ignoreNormCorrectionFlag = false;
  m_ignoreRandCorrectionFlag = false;
  m_ignoreScatCorrectionFlag = false;
  m_ignoreDecaCorrectionFlag = false;
  m_ignoreBratCorrectionFlag = false;
  m_ignoreFdurCorrectionFlag = false;
  m_ignoreCaliCorrectionFlag = false;
  m_nbRespGates = 1;
  m_nbRespBasisFunctions = -1;
  m2p_respBasisFunctions = NULL;
  m_respBasisFunctionsFile = "";
  m_respStaticFlag = false;
  m_nbCardGates = 1;
  m_nbCardBasisFunctions = -1;
  m2p_cardBasisFunctions = NULL;
  m_cardBasisFunctionsFile = "";
  m_cardStaticFlag = false;
  m_nbThreadsForProjection = 1;
  m_nbThreadsForImageComputation = 1;
  m_mpiRank = 0;
  m_mpiSize = 1;
  m_nbBeds = -1;
  mp_bedPositions = NULL;
  m_providedBedPosition = false;
  m_verbose = 0;
  m_checked = false;
  m_initialized = false;
  m_dynRecoTypeFlag = STATIC_RECO;
  // Allocate Dynamic data manager object
  mp_DynamicDataManager = new oDynamicDataManager();
  m_nbMultiModalImages = 0;
  m_lambda = -1.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oImageDimensionsAndQuantification::~oImageDimensionsAndQuantification()
{
  // Free frame duration
  if (m2p_frameDurationsInMs)
  {
    for (int bed=0; bed<m_nbBeds; bed++) if (m2p_frameDurationsInMs[bed]) delete m2p_frameDurationsInMs[bed];
    delete m2p_frameDurationsInMs;
  }
  // Free frame time start
  if (m2p_frameTimeStartInMs)
  {
    for (int bed=0; bed<m_nbBeds; bed++) if (m2p_frameTimeStartInMs[bed]) delete m2p_frameTimeStartInMs[bed];
    delete m2p_frameTimeStartInMs;
  }
  // Free frame time stop
  if (m2p_frameTimeStopInMs)
  {
    for (int bed=0; bed<m_nbBeds; bed++) if (m2p_frameTimeStopInMs[bed]) delete m2p_frameTimeStopInMs[bed];
    delete m2p_frameTimeStopInMs;
  }
  // Free quantification factors
  if (m3p_quantificationFactors)
  {
    for (int bed=0; bed<m_nbBeds; bed++)
    {
      if (m3p_quantificationFactors[bed])
      {
        for (int fr=0; fr<m_nbTimeFrames; fr++)
          if (m3p_quantificationFactors[bed][fr]) delete m3p_quantificationFactors[bed][fr];
      }
      delete m3p_quantificationFactors[bed];
    }
    delete m3p_quantificationFactors;
  }
  // Free the bed positions
  if (mp_bedPositions) free(mp_bedPositions);
  // Delete the dynamic data manager
  delete mp_DynamicDataManager;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::SetFlipOut(const string& a_flipOut)
{
  // If empty, then put all flags to false and return
  if (a_flipOut=="")
  {
    m_flipOutX = false;
    m_flipOutY = false;
    m_flipOutZ = false;
    return 0;
  }
  // Test all possible settings !
  if (a_flipOut=="x" || a_flipOut=="X")
  {
    m_flipOutX = true;
    m_flipOutY = false;
    m_flipOutZ = false;
  }
  else if (a_flipOut=="y" || a_flipOut=="Y")
  {
    m_flipOutX = false;
    m_flipOutY = true;
    m_flipOutZ = false;
  }
  else if (a_flipOut=="z" || a_flipOut=="Z")
  {
    m_flipOutX = false;
    m_flipOutY = false;
    m_flipOutZ = true;
  }
  else if (a_flipOut=="xy" || a_flipOut=="yx" || a_flipOut=="XY" || a_flipOut=="YX")
  {
    m_flipOutX = true;
    m_flipOutY = true;
    m_flipOutZ = false;
  }
  else if (a_flipOut=="zy" || a_flipOut=="yz" || a_flipOut=="ZY" || a_flipOut=="YZ")
  {
    m_flipOutX = false;
    m_flipOutY = true;
    m_flipOutZ = true;
  }
  else if (a_flipOut=="xz" || a_flipOut=="zx" || a_flipOut=="XZ" || a_flipOut=="ZX")
  {
    m_flipOutX = true;
    m_flipOutY = false;
    m_flipOutZ = true;
  }
  else if ( a_flipOut=="xyz" || a_flipOut=="xzy" || a_flipOut=="yxz" || a_flipOut=="yzx" || a_flipOut=="zxy" || a_flipOut=="zyx" ||
            a_flipOut=="XYZ" || a_flipOut=="XZY" || a_flipOut=="YXZ" || a_flipOut=="YZX" || a_flipOut=="ZXY" || a_flipOut=="ZYX" )
  {
    m_flipOutX = true;
    m_flipOutY = true;
    m_flipOutZ = true;
  }
  // Otherwise, throw an error
  else
  {
    Cerr("***** oImageDimensionsAndQuantification::SetFlipOut() -> Output flip settings is incorrect !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::SetNbThreads(const string& a_nbThreads)
{
  // The number of threads can be given as two separated numbers to distinguish between the number of threads for projection computation
  // and the number of threads for image computation. Be careful that the thread-safe images are allocated with respect to the number
  // of projection threads; the number of threads used for image computation are strictly used for computation over images (i.e. when
  // using a loop over voxels.

  // Look for a comma and check that there is only one comma without missing parameters
  size_t first_comma = a_nbThreads.find_first_of(",");
  size_t last_comma = a_nbThreads.find_last_of(",");
  if (first_comma!=last_comma || first_comma==0 || first_comma==a_nbThreads.length()-1)
  {
    Cerr("***** oImageDimensionsAndQuantification::SetNbThreads() -> Wrong syntax in the thread parameters ! See help." << endl);
    return 1;
  }

  // Case for a single parameter
  if (first_comma==string::npos)
  {
    m_nbThreadsForProjection = atoi( a_nbThreads.c_str() );
    m_nbThreadsForImageComputation = m_nbThreadsForProjection;
  }
  // Case for two parameters
  else
  {
    m_nbThreadsForProjection = atoi( (a_nbThreads.substr(0,first_comma)).c_str() );
    m_nbThreadsForImageComputation = atoi( (a_nbThreads.substr(first_comma+1)).c_str() );
  }

  // Checks for negative threads numbers
  if (m_nbThreadsForProjection<0)
  {
    Cerr("***** oImageDimensionsAndQuantification::SetNbThreads() -> Negative number of threads provided for projection computation !" << endl);
    return 1;
  }
  if (m_nbThreadsForImageComputation<0)
  {
    Cerr("***** oImageDimensionsAndQuantification::SetNbThreads() -> Negative number of threads provided for image computation !" << endl);
    return 1;
  }

  #ifdef CASTOR_OMP
  // If number of threads is 0, we automatically determine the maximum number of threads using OMP functions
  if (m_nbThreadsForProjection==0) m_nbThreadsForProjection = omp_get_max_threads();
  if (m_nbThreadsForImageComputation==0) m_nbThreadsForImageComputation = omp_get_max_threads();
  // At this step, set by default the number of threads for image operations
  omp_set_num_threads(m_nbThreadsForImageComputation);
  #endif

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageDimensionsAndQuantification::SetDefault()
{
  // Number of threads
  m_nbThreadsForProjection = 1;
  m_nbThreadsForImageComputation = 1;
  // MPI stuff
  m_mpiSize = 1;
  m_mpiRank = 0;
  // Unlock
  m_checked = true;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::CheckParameters()
{
  // Number of threads
  if (m_nbThreadsForProjection<=0)
  {
    Cerr("***** oImageDimensionsAndQuantification::CheckParameters() -> Should provide a strictly positive number of threads !" << endl);
    return 1;
  }
  // TODO: authorize -fov and -vox without -dim, checking that the number of voxels obtained is an integer
  // Number of voxels
  if (m_nbVoxX<=0 || m_nbVoxY<=0 || m_nbVoxZ<=0)
  {
    Cerr("***** oImageDimensionsAndQuantification::CheckParameters() -> Should provide strictly positive number of voxels !" << endl);
    return 1;
  }
  // When FOV and voxel sizes are both not provided
  if ((m_voxSizeX<=0. || m_voxSizeY<=0. || m_voxSizeZ<=0.) && (m_fovSizeX<=0. || m_fovSizeY<=0. || m_fovSizeZ<=0.))
  {
    Cerr("***** oImageDimensionsAndQuantification::CheckParameters() -> Should provide strictly positive voxel or FOV dimensions !" << endl);
    return 1;
  }
  // When both FOV and voxel sizes are provided
  if (m_voxSizeX>0. && m_voxSizeY>0. && m_voxSizeZ>0. && m_fovSizeX>0. && m_fovSizeY>0. && m_fovSizeZ>0.)
  {
    Cerr("***** oImageDimensionsAndQuantification::CheckParameters() -> Both FOV and voxels dimensions provided, should not provide both !" << endl);
    return 1;
  }
  // Check output transaxial FOV which is in percent (0 means we do not apply any FOV masking)
  if (m_fovOutPercent<0.)
  {
    Cerr("***** oImageDimensionsAndQuantification::CheckParameters() -> Output transaxial FOV percentage must be strictly positive !" << endl);
    return 1;
  }
  // Check that the number of axial slices to be masked is between 0 and dimZ/2
  if (m_nbSliceOutMask<0 || m_nbSliceOutMask>m_nbVoxZ/2)
  {
    Cerr("***** oImageDimensionsAndQuantification::CheckParameters() -> Number of output axial slices to be masked is incorrectly set !" << endl);
    return 1;
  }
  // Check whether the DynamicDataManager object has been instanciated or not
  if (!mp_DynamicDataManager)
  {
    Cerr("***** oImageDimensionsAndQuantification::CheckParameters() -> Error : DynamicDataManager object not initialized !" << endl);
    return 1;
  }
  // Unlock
  m_checked = true;
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::Initialize()
{
  // Mandatory check
  if (!m_checked)
  {
    Cerr("***** oImageDimensionsAndQuantification::Initialize() -> Cannot initialize before a call to CheckParameters() !" << endl);
    return 1;
  }

  // Verbose
  if (m_verbose>=2) Cout("oImageDimensionsAndQuantification::Initialize() -> Initialize image dimensions, basis functions and quantification" << endl);

  // Precompute number of voxels in a slice (XY) and total number of voxels (XYZ)
  m_nbVoxXY = m_nbVoxX * m_nbVoxY;
  m_nbVoxXYZ = m_nbVoxXY * m_nbVoxZ;

  // If FOV dimensions are provided, then compute the voxel dimensions
  if (m_fovSizeX>0. && m_fovSizeY>0. && m_fovSizeZ>0.)
  {
    m_voxSizeX = m_fovSizeX / ((FLTNB)m_nbVoxX);
    m_voxSizeY = m_fovSizeY / ((FLTNB)m_nbVoxY);
    m_voxSizeZ = m_fovSizeZ / ((FLTNB)m_nbVoxZ);
  }
  // If voxel dimensions are provided, then compute the voxel dimensions
  else if (m_voxSizeX>0. && m_voxSizeY>0. && m_voxSizeZ>0.)
  {
    m_fovSizeX = m_voxSizeX * ((FLTNB)m_nbVoxX);
    m_fovSizeY = m_voxSizeY * ((FLTNB)m_nbVoxY);
    m_fovSizeZ = m_voxSizeZ * ((FLTNB)m_nbVoxZ);
  }

  // Deal with frame list and quantification factors
  if (InitializeFramingAndQuantification())
  {
    Cerr("***** oImageDimensionsAndQuantification::Initialize() -> A problem occurred while initializing framing and quantification tabs !" << endl);
    return 1;
  }

  // Deal with time basis functions
  /*
  if (InitializeTimeBasisFunctions())
  {
    Cerr("***** oImageDimensionsAndQuantification::Initialize() -> A problem occurred while initializing time basis functions !" << endl);
    return 1;
  }

  // Deal with respiratory basis functions
  if (InitializeRespBasisFunctions())
  {
    Cerr("***** oImageDimensionsAndQuantification::Initialize() -> A problem occurred while initializing respiratory basis functions !" << endl);
    return 1;
  }

  // Deal with cardiac basis functions
  if (InitializeCardBasisFunctions())
  {
    Cerr("***** oImageDimensionsAndQuantification::Initialize() -> A problem occurred while initializing cardiac basis functions !" << endl);
    return 1;
  }
//*/
  // Deal with ignored corrections
  if (InitializeIgnoredCorrections())
  {
    Cerr("***** oImageDimensionsAndQuantification::Initialize() -> A problem occurred while initializing ignored corrections !" << endl);
    return 1;
  }

  // Verbose
  if (m_verbose>=2)
  {
    Cout("  --> Image dimensions: [" << m_nbVoxX << ";" << m_nbVoxY << ";" << m_nbVoxZ << "] voxels of [" << m_voxSizeX << ";" << m_voxSizeY << ";" << m_voxSizeZ << "] mm3" << endl);
    Cout("  --> FOV size: [" << m_fovSizeX << ";" << m_fovSizeY << ";" << m_fovSizeZ << "] mm3" << endl);
    if (m_nbTimeFrames>1)
    {
      if (m_timeStaticFlag) Cout("  --> Time frames: " << m_nbTimeFrames << endl);
      else Cout("  --> Time frames: " << m_nbTimeFrames << " | Time basis functions: " << m_nbTimeBasisFunctions << endl);
    }
    if (m_respBasisFunctionsFile!="")
    {
      if (m_respStaticFlag) Cout("  --> Respiratory gates: " << m_nbRespGates << endl);
      else Cout("  --> Respiratory gates: " << m_nbRespGates << " | Respiratory basis functions: " << m_nbRespBasisFunctions << endl);
    }
    if (m_cardBasisFunctionsFile!="")
    {
      if (m_cardStaticFlag) Cout("  --> Cardiac gates: " << m_nbCardGates << endl);
      else Cout("  --> Cardiac gates: " << m_nbCardGates << " | Cardiac basis functions: " << m_nbCardBasisFunctions << endl);
    }
    if (m_nbThreadsForProjection>1 || m_nbThreadsForImageComputation>1)
    {
      if (m_nbThreadsForImageComputation==m_nbThreadsForProjection) Cout("  --> Number of parallel threads: " << m_nbThreadsForProjection << endl);
      else Cout("  --> Number of parallel threads for projection / image computation: [" << m_nbThreadsForProjection << "/" << m_nbThreadsForImageComputation << "]" << endl);
    }
  }

  // Initialized
  m_initialized = true;

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageDimensionsAndQuantification::CheckNumberOfProjectionThreadsConsistencyWithDataFileSize(vDataFile** a2p_DataFile)
{
  // Loop on all datafiles
  for (int bed=0; bed<m_nbBeds; bed++)
  {
    // If the number of events in this datafile is below the number of projection threads, then reduce the number of projection threads
    // to this number
    if (a2p_DataFile[bed]->GetSize()<m_nbThreadsForProjection)
    {
      m_nbThreadsForProjection = a2p_DataFile[bed]->GetSize();
      Cerr("!!!!! oImageDimensionsAndQuantification::CheckNumberOfProjectionThreadsConsistencyWithDataFileSize() !!!!!" << endl);
      Cerr("      --> The number of projection threads was reduced to the provided datafile's number of events: " << m_nbThreadsForProjection << endl);
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::DealWithBedPositions(vDataFile** a2p_DataFile)
{
  // Note that the consistency between bed positions should have already been tested when this function is called.

  // First particular case for a single bed
  if (m_nbBeds==1)
  {
    // Here we just copy the relative bed position if provided by the datafile header. It will be ignored in the projection when
    // there is only one bed position, so no worries.
    mp_bedPositions[0] = a2p_DataFile[0]->GetRelativeBedPosition();
    // Flag to say that the bed positions were provided from the datafile
    m_providedBedPosition = a2p_DataFile[0]->GetBedPositionFlag();
    // We directly exit here, because below in this function, if the relative bed position has not been provided in the datafile header,
    // then the scanner default bed displacement will be used, however, this value is not mandatory as some scanners may not be designed
    // for wholebody.
    return 0;
  }

  // Case where the bed relative positions were provided for all the datafiles
  if (a2p_DataFile[0]->GetBedPositionFlag())
  {
    // In this case we have to compute the center of mass of the provided bed positions, and move it to the 0-center of the Z axis
    FLTNB center = 0.;
    // Compute the center of mass
    for (int bed=0; bed<m_nbBeds; bed++) center += a2p_DataFile[bed]->GetRelativeBedPosition();
    center /= ((FLTNB)m_nbBeds);
    // Compute the new bed positions from their relative values to recenter the mass at 0
    for (int bed=0; bed<m_nbBeds; bed++) mp_bedPositions[bed] = a2p_DataFile[bed]->GetRelativeBedPosition() - center;
    // Flag to say that the bed positions were provided from the datafile
    m_providedBedPosition = true;
  }

  // Case where the scanner default bed displacement is used
  else
  {
    // Get the scanner object in use
    vScanner* ap_Scanner = sScannerManager::GetInstance()->GetScannerObject();
    // Check that the default bed displacement is more than 0. in the scanner
    if (ap_Scanner->GetDefaultBedDisplacementInMm()<=0.)
    {
      Cerr("***** oImageDimensionsAndQuantification::DealWithBedPositions() -> Bed displacement between two successive bed positions must be strictly positive !" << endl);
      return 1;
    }
    // Loop on the bed positions
    for (int bed=0; bed<m_nbBeds; bed++)
    {
      // Compute bed offset (remember here that the 0. on the Z axis is at the center of the whole image)
      FLTNB bed_offset = 0.;
      // For odd numbers of bed positions
      if (m_nbBeds%2==1) bed_offset = ((FLTNB)( bed-m_nbBeds/2 )) * ap_Scanner->GetDefaultBedDisplacementInMm();
      // For even numbers of bed positions
      else bed_offset = (((FLTNB)( bed-m_nbBeds/2 )) + 0.5) * ap_Scanner->GetDefaultBedDisplacementInMm();
      // Record the value
      mp_bedPositions[bed] = bed_offset;
    }
    // Flag to say that the bed positions were not provided from the datafile
    m_providedBedPosition = false;
  }

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("oImageDimensionsAndQuantification::DealWithBedPositions() -> Use following relative bed positions:" << endl);
    for (int bed=0; bed<m_nbBeds; bed++) Cout("  --> Bed " << bed << " | Relative axial position: " << mp_bedPositions[bed] << " mm" << endl);
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::InitializeFramingAndQuantification()
{
// TODO: for whole body 4D dynamic acquisitions, will need to read a set of frames for each bed where the number of frames and frame durations
//       should be the same for all beds, except which frame is dead and time start and stop also... Have to think about and manage this.
  // This function requires that the number of beds be set
  if (m_nbBeds<=0)
  {
    Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Number of beds must be set before setting the framing of the acquisition !" << endl);
    return 1;
  }

  // --------------------------------------------------
  // Allocate time start, stop and duration, to 1 frame
  m2p_frameDurationsInMs    = (uint32_t**)malloc(m_nbBeds*sizeof(uint32_t*));
  m2p_frameTimeStartInMs    = (uint32_t**)malloc(m_nbBeds*sizeof(uint32_t*));
  m2p_frameTimeStopInMs     = (uint32_t**)malloc(m_nbBeds*sizeof(uint32_t*));
  m3p_quantificationFactors = (FLTNB***)malloc(m_nbBeds*sizeof(FLTNB**));
  for (int bed=0; bed<m_nbBeds; bed++)
  {
    m2p_frameDurationsInMs[bed]    = (uint32_t*)malloc(1*sizeof(uint32_t));
    m2p_frameTimeStartInMs[bed]    = (uint32_t*)malloc(1*sizeof(uint32_t));
    m2p_frameTimeStopInMs[bed]     = (uint32_t*)malloc(1*sizeof(uint32_t));
  }

  // ---------------------------------------
  // Particular case where the list is empty
  if (m_frameList=="")
  {
    // Set number of frames to 1
    m_nbTimeFrames = 1;
    // Init them to 0 (the vDataFile will set them to the appropriate value later through the oImageDimensionsAndQuantification::SetAcquisitionTime() function)
    for (int bed=0; bed<m_nbBeds; bed++)
    {
      m2p_frameDurationsInMs[bed][0] = 0;
      m2p_frameTimeStartInMs[bed][0] = 0;
      m2p_frameTimeStopInMs[bed][0]  = 0;
      // Allocate quantification factors and set them to 1.
      m3p_quantificationFactors[bed] = (FLTNB**)malloc(1*sizeof(FLTNB*));
      m3p_quantificationFactors[bed][0] = (FLTNB*)malloc(m_nbRespGates*m_nbCardGates*sizeof(FLTNB));
      for (int g=0; g<m_nbRespGates*m_nbCardGates; g++) m3p_quantificationFactors[bed][0][g] = 1.;
    }
    // Exit the function
    return 0;
  }

  // --------------------------------------------------------------------
  // Otherwise, get the parameter as a non-constant string and process it
  string frame_list = m_frameList;
  
  // Declare frame start and duration
  HPFLTNB frame_start ;
  HPFLTNB frame_duration ;
  // Units are defaulted to seconds, unless minutes are specified within the frame decleration.
  bool frame_start_inMinutes = false;
  bool frame_duration_inMinutes = false;
  // Init number of frames
  m_nbTimeFrames = 0;
  // Loop on all commas and colons found in the list
  size_t comma_pos = 0;
  size_t colon_pos = 0;
  size_t unit_pos = 0;
  while ((comma_pos=frame_list.find_first_of(","))!=string::npos)
  {
    // Set default values to 0
    frame_start = (HPFLTNB) 0.;
    frame_duration = (HPFLTNB) 0.;
    // By default set/reset units to seconds
    frame_start_inMinutes = false;
    frame_duration_inMinutes = false;
    // Increment number of frames by 1
    m_nbTimeFrames++;
    // Realloc the time start, stop and duration arrays for all beds
    for (int bed=0; bed<m_nbBeds; bed++)
    {
      m2p_frameDurationsInMs[bed] = (uint32_t*)realloc( m2p_frameDurationsInMs[bed], m_nbTimeFrames*sizeof(uint32_t) );
      m2p_frameTimeStartInMs[bed] = (uint32_t*)realloc( m2p_frameTimeStartInMs[bed], m_nbTimeFrames*sizeof(uint32_t) );
      m2p_frameTimeStopInMs[bed]  = (uint32_t*)realloc( m2p_frameTimeStopInMs[bed] , m_nbTimeFrames*sizeof(uint32_t) );
    }
    // Extract the parameter ("phrase") before the first comma
    string param = frame_list.substr(0,comma_pos);
    // Check if the parameter is empty.
    if (param.empty())
    {
      Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Null framing definition detected !" << endl);
      return 1;
    }
    // Search for a colon within the param (this means that time 'start:duration' are both provided)
    if ((colon_pos = param.find_first_of(":")) != string::npos)
    {
      // Getting frame start declaration
      string param_start = param.substr(0, colon_pos);
      // Search for definition of units within the frame start declaration
      if ((unit_pos = param_start.find("s"))!= string::npos)
      {
        // If seconds unit found - remove from param_start
        param_start.erase(unit_pos);
      }
      else if ((unit_pos = param_start.find("m"))!= string::npos)
      {
        // If minutes unit found - remove from param_start
        param_start.erase(unit_pos);
        frame_start_inMinutes = true;
      }
      // Checking that frame start declaration is not empty
      if (param_start.empty())
      {
        Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Null framing definition detected !" << endl);
        return 1;
      }
      // Getting frame duration declaration
      string param_duration = param.substr(colon_pos + 1, comma_pos);
      // Search for definition of units
      if ((unit_pos = param_duration.find("s"))!= string::npos)
      {
        param_duration.erase(unit_pos);
      }
      else if ((unit_pos = param_duration.find("m"))!= string::npos)
      {
        param_duration.erase(unit_pos);
        frame_duration_inMinutes = true;
      }
      // Checking that frame duration declaration is not empty
      if (param_duration.empty())
      {
        Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Null framing duration detected !" << endl);
        return 1;
      }
      // Assign numerical values to frame start and duration
      frame_start = (HPFLTNB) atof(param_start.c_str());
      frame_duration = (HPFLTNB) atof(param_duration.c_str());
      // If units declared in minutes - convert to seconds
      if (frame_start_inMinutes) frame_start *= (HPFLTNB)60.;
      if (frame_duration_inMinutes) frame_duration *= (HPFLTNB)60.;
      // Check for negative or null frame duration
      if (frame_duration<=0.)
      {
        Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Negative or null frame duration detected !" << endl);
        return 1;
      }
      // Check for negative or null frame start time
      if (frame_start<0.)
      {
        Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Negative or null frame start detected !" << endl);
        return 1;
      }
      // Affect time start and duration; Conversion to milliseconds and cast to unit32_t.
      for (int bed=0; bed<m_nbBeds; bed++)
      {
        m2p_frameDurationsInMs[bed][m_nbTimeFrames-1] = ((uint32_t)(frame_duration*1000.));
        m2p_frameTimeStartInMs[bed][m_nbTimeFrames-1] = ((uint32_t)(frame_start*1000.));
      }
    }
    // Case where frame duration has not been declared (no colon found)
    else
    {
      // Search for definition of units within the frame start declaration
      if ((unit_pos = param.find("s"))!= string::npos)
      {
        // If seconds unit found - remove from param_start
        param.erase(unit_pos);
      }
      else if ((unit_pos = param.find("m"))!= string::npos)
      {
        // If minutes unit found - remove from param_start
        param.erase(unit_pos);
        frame_start_inMinutes=true;
      }
      // Check that frame start is not empty
      if (param.empty())
      {
        Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Null framing definition detected !" << endl);
        return 1;
      }
      // Assign value
      frame_start = (HPFLTNB) atof(param.c_str());
      // If units declared in minutes - convert to seconds
      if (frame_start_inMinutes) frame_start *= (HPFLTNB)60.;
      // Check for negative or null frame start time
      if (frame_start<0.)
      {
        Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Negative or null frame start detected !" << endl);
        return 1;
      }
      // Affect time start
      for (int bed=0; bed<m_nbBeds; bed++)
      {
        // Frame duration not provided, we set 0 temporarily; it will be computed after reading all frames, and set to fit until the next frame start
        m2p_frameDurationsInMs[bed][m_nbTimeFrames-1] = 0;
        m2p_frameTimeStartInMs[bed][m_nbTimeFrames-1] = ((uint32_t)(frame_start*1000.));
      }
    }
    // Remove the evaluated parameter "phrase"
    frame_list = frame_list.substr(comma_pos+1);
  }

  // --- Last parameter to evaluate and extract ---
  // Last frame must always have a declaration of frame duration

  // Set default values to 0
  frame_start = (HPFLTNB) 0.;
  frame_duration = (HPFLTNB) 0.;
  // By default units set to seconds
  frame_start_inMinutes = false;
  frame_duration_inMinutes = false;
  // Case with the duration (here it is mandatory for the last parameter)
  string param = frame_list;
  if ((colon_pos = param.find_first_of(":")) != string::npos)
  {
    // Increment number of frames
    m_nbTimeFrames++;
    // Realloc the time start, stop and duration
    for (int bed=0; bed<m_nbBeds; bed++)
    {
      m2p_frameDurationsInMs[bed]    = (uint32_t*)realloc(m2p_frameDurationsInMs[bed],m_nbTimeFrames*sizeof(uint32_t));
      m2p_frameTimeStartInMs[bed]    = (uint32_t*)realloc(m2p_frameTimeStartInMs[bed],m_nbTimeFrames*sizeof(uint32_t));
      m2p_frameTimeStopInMs[bed]     = (uint32_t*)realloc(m2p_frameTimeStopInMs[bed] ,m_nbTimeFrames*sizeof(uint32_t));
    }
    // Getting frame start declaration
    string param_start = param.substr(0, colon_pos);
    // Search for definition of units within the frame start declaration
    if ((unit_pos = param_start.find("s"))!= string::npos)
    {
      // if minutes unit found - remove from param_start
      param_start.erase(unit_pos);
    }
    else if ((unit_pos = param_start.find("m"))!= string::npos)
    {
      // if seconds unit found - remove from param_start
      param_start.erase(unit_pos);
      frame_start_inMinutes = true;
    }
    // Check that start frame is not empty
    if (param_start.empty())
    {
      Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Null framing definition detected !" << endl);
      return 1;
    }
    // Getting frame duration declaration
    string param_duration = param.substr(colon_pos + 1, param.size());
    // Search for definition of units
    if ((unit_pos = param_duration.find("s"))!= string::npos)
    {
      param_duration.erase(unit_pos);
    }
    else if ((unit_pos = param_duration.find("m"))!= string::npos)
    {
      param_duration.erase(unit_pos);
      frame_duration_inMinutes=true;
    }
    // Checking that frame duration is not empty
    if (param_duration.empty())
    {
      Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Null framing duration detected !" << endl);
      return 1;
    }
    // Assign numerical values to frame start and duration
    frame_start = (HPFLTNB) atof(param_start.c_str());
    frame_duration = (HPFLTNB) atof(param_duration.c_str());
    // If units declared in minutes - convert to seconds
    if (frame_start_inMinutes) frame_start *= (HPFLTNB)60.;
    if (frame_duration_inMinutes) frame_duration *= (HPFLTNB)60.;
    // Check for negative or null duration
    if (frame_duration<=0.)
    {
      Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Found a duration which is negative or null !" << endl);
      return 1;
    }
    // Check for negative frame start
    if (frame_start<0.)
    {
      Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Negative or null frame start detected !" << endl);
      return 1;
    }
    // Affect time start, stop and duration
    for (int bed=0; bed<m_nbBeds; bed++)
    {
      m2p_frameDurationsInMs[bed][m_nbTimeFrames-1] = ((uint32_t)(frame_duration*1000.));
      m2p_frameTimeStartInMs[bed][m_nbTimeFrames-1] = ((uint32_t)(frame_start*1000.));
    }

  }
  // Otherwise, if no colon was found in last frame, give an error as the last frame must always have declaration of duration
  else
  {
    Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Last frame duration has not been provided !" << endl);
    return 1;
  }

  // -----------------------------------------------
  // Compute the frame durations for the frames it has not been declared, skipping the last frame as its duration is always declared
  for (int frm=0; frm<m_nbTimeFrames-1; frm++)
  {
    // Duration was set to 0 when not provided
    if (m2p_frameDurationsInMs[0][frm] == 0)
    {
      for (int bed=0; bed<m_nbBeds; bed++)
      {
        m2p_frameDurationsInMs[bed][frm] = m2p_frameTimeStartInMs[bed][frm+1] - m2p_frameTimeStartInMs[bed][frm];
      }
    }
  }

  // -----------------------------------------------
  // Allocate the frame end times for all frames
  for (int frm=0; frm<m_nbTimeFrames; frm++)
  {
    for (int bed=0; bed<m_nbBeds; bed++)
    {
      m2p_frameTimeStopInMs[bed][frm] = m2p_frameTimeStartInMs[bed][frm] + m2p_frameDurationsInMs[bed][frm];
    }
  }

  // -------------------------------------------------
  // Checking for Overlap of frames
  // Checking that for every single frame, the start time and end time do not fall within any other frame
  // TODO: checking a single bed at the moment - as definition of times for all beds is the same

  for (int frmch=0; frmch<m_nbTimeFrames; frmch++)
  {
    for (int frm=frmch+1; frm<m_nbTimeFrames; frm++)
    {
      if (m2p_frameTimeStartInMs[0][frm] < m2p_frameTimeStopInMs[0][frmch])
      {
        Cerr("***** oImageDimensionsAndQuantification::InitializeFramingAndQuantification() -> Illegal frame overlap detected between frames: " << frmch+1 << " and "<< frm+1 << endl);
        return 1;
      }
    }
  }

  // ----------------------------------------------
  // Allocate and affect the quantification factors
  for (int bed=0; bed<m_nbBeds; bed++)
  {
    m3p_quantificationFactors[bed] = (FLTNB**)malloc(m_nbTimeFrames*sizeof(FLTNB*));
    for (int fr=0; fr<m_nbTimeFrames; fr++)
    {
      m3p_quantificationFactors[bed][fr] = (FLTNB*)malloc(m_nbRespGates*m_nbCardGates*sizeof(FLTNB));
      for (int g=0; g<m_nbRespGates*m_nbCardGates; g++) m3p_quantificationFactors[bed][fr][g] = 1.;
    }
  }

  // ---
  // End
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::InitializeIgnoredCorrections()
{
  // If the m_ignoredCorrectionsList is empty, force all booleans to false and leave this function
  if (m_ignoredCorrectionsList=="")
  {
    m_ignoreAttnCorrectionFlag = false;
    m_ignoreNormCorrectionFlag = false;
    m_ignoreRandCorrectionFlag = false;
    m_ignoreScatCorrectionFlag = false;
    m_ignoreDecaCorrectionFlag = false;
    m_ignoreBratCorrectionFlag = false;
    m_ignoreFdurCorrectionFlag = false;
    m_ignoreCaliCorrectionFlag = false;
    return 0;
  }

  // Then, count the number of commas in the m_ignoredCorrectionsList to know the number of keywords
  size_t nb_keywords = count(m_ignoredCorrectionsList.begin(), m_ignoredCorrectionsList.end(), ',') + 1;

  // Read the keywords
  string *p_keywords = new string[nb_keywords];

  if (ReadStringOption(m_ignoredCorrectionsList, p_keywords, nb_keywords, ",", "Ignored corrections list"))
  {
    Cerr("***** oImageDimensionsAndQuantification::InitializeIgnoredCorrections() -> An error occurred while reading the list of ignored corrections !" << endl);
    return 1;
  }

  // Process them
  for (size_t k=0; k<nb_keywords; k++)
  {
    // Test the different keywords
    if (p_keywords[k]=="attn") m_ignoreAttnCorrectionFlag = true;
    else if (p_keywords[k]=="norm") m_ignoreNormCorrectionFlag = true;
    else if (p_keywords[k]=="scat") m_ignoreScatCorrectionFlag = true;
    else if (p_keywords[k]=="rand") m_ignoreRandCorrectionFlag = true;
    else if (p_keywords[k]=="deca") m_ignoreDecaCorrectionFlag = true;
    else if (p_keywords[k]=="brat") m_ignoreBratCorrectionFlag = true;
    else if (p_keywords[k]=="fdur") m_ignoreFdurCorrectionFlag = true;
    else if (p_keywords[k]=="cali") m_ignoreCaliCorrectionFlag = true;
    else
    {
      Cerr("***** oImageDimensionsAndQuantification::InitializeIgnoredCorrections() -> Unknown keyword '" << p_keywords[k] << "' in the provided ignored corrections list !" << endl);
      return 1;
    }
  }

  delete[] p_keywords;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::SetAcquisitionTime(int a_bed, FLTNB a_timeStartInSec, FLTNB a_durationInSec, string a_gateListDurationInSec)
{
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageDimensionsAndQuantification::SetAcquisitionTime() -> Object not initialized !" << endl);
    return 1;
  }
  
  // TODO TM : Once the dynamic format will be updated (one file by frame/gate)
  // The gate duration will have to be read from the datafile header, and the quantification factors must be updated here
  // (Yet to see the most practical way to manage information with several datafiles
  // Hete, the quantification factors for all frames and gates are simultaneously updated)
  //
  // For now, just one factor by gate is considered, for data without dynamic frame (must be frame dependent /!!!!!!\)
  // Parse the list containing the gate durations 
  // Local array to recover gate parameters
  FLTNB* pgate_duration_sec = new FLTNB[ m_nbRespGates*m_nbCardGates ];
  
  // First check that gate duration have not been provided in BOTH the datafile header and the gating configuration file, the later being the current mandatory way (make your mind!)
  if(a_gateListDurationInSec != "" && mp_DynamicDataManager->GateDurationProvided())
  {
    Cerr("***** oImageDimensionsAndQuantification::SetAcquisitionTime() -> Gate durations have been initialized in both the datafile header and the gating configuration file!" << endl);
    Cerr("*****                                                            Only one must be used (preferably the gating configuration file) " << endl);
    return 1;  
  }
  
  if(a_gateListDurationInSec != "")
  {
    if (ReadStringOption(a_gateListDurationInSec,
                         pgate_duration_sec,
                         m_nbRespGates*m_nbCardGates,
                         ",",
                         "Gate duration (s)"))
    {
      Cerr("***** oImageDimensionsAndQuantification::SetAcquisitionTime() -> Failed to correctly read the following list of gate durations (datafile header): !" << endl);
      Cerr("*****                                                            "<<a_gateListDurationInSec << endl);
      Cerr("*****                                                            "<<m_nbRespGates*m_nbCardGates<<" parameters were expected (1 for each gate)" << endl);
      return 1;
    }
  }
  // No gate, initialize with the acquisition duration
  else
  {
    // Initialize with acquisition duration, unless gate duration have been provided and quantification factors already integrate time quantification correction 
    // TODO TM: This has to be deleted when new dynamic datafile format will be implemented (gate duration directly in datafile header)
    for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
      pgate_duration_sec[ g ]= mp_DynamicDataManager->GateDurationProvided() ? 1 : a_durationInSec;
  }
  
  // If we have only one frame and that, timeStart, timeStop and duration are 0, then it means that nothing has been specified yet.
  // (this was done in the InitializeFramingAndQuantification() function)
  if (m_nbTimeFrames==1 && m2p_frameDurationsInMs[a_bed][0]==0 && m2p_frameTimeStartInMs[a_bed][0]==0 && m2p_frameTimeStopInMs[a_bed][0]==0)
  {
    m2p_frameDurationsInMs[a_bed][0] = ((uint32_t)(a_durationInSec*1000.));
    m2p_frameTimeStartInMs[a_bed][0] = ((uint32_t)(a_timeStartInSec*1000.));
    m2p_frameTimeStopInMs[a_bed][0]  = m2p_frameTimeStartInMs[a_bed][0] + m2p_frameDurationsInMs[a_bed][0];
    // Apply frame's duration onto quantification factors
    if (!m_ignoreFdurCorrectionFlag)
    {
      for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
        m3p_quantificationFactors[a_bed][0][g] /= pgate_duration_sec[ g ];
    }
  }
  // Otherwise it was already specified, so we just update the quantification factors and exit
  else
  {
    if (!m_ignoreFdurCorrectionFlag 
     && !mp_DynamicDataManager->GateDurationProvided() ) // Durations provided in gate configuration file. No need to update quantification factors.
    {
      for (int fr=0; fr<m_nbTimeFrames; fr++) for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
        m3p_quantificationFactors[a_bed][fr][g] /= (((FLTNB)(m2p_frameDurationsInMs[a_bed][fr]))/1000.);
    }
  }
  
  // Verbose
  if (m_verbose>=2 && (a_bed==m_nbBeds-1))
  {
    // Case 1: a static single bed
    if (m_nbTimeFrames==1 && m_nbBeds==1)
    {
      Cout("oImageDimensionsAndQuantification::SetAcquisitionTime() -> Static single bed acquisition with duration [ " << GetFrameTimeStartInSec(0,0) << " : "
           << GetFrameTimeStopInSec(0,0) << " ] seconds" << endl);
    }
    // Case 2: a static multi beds
    else if (m_nbTimeFrames==1 && m_nbBeds>1)
    {
      Cout("oImageDimensionsAndQuantification::SetAcquisitionTime() -> Static " << m_nbBeds << " beds acquisition with following bed durations:" << endl);
      for (int bed=0; bed<m_nbBeds; bed++)
        Cout("  --> Bed " << bed+1 << " with duration [ " << GetFrameTimeStartInSec(bed,0) << " : " << GetFrameTimeStopInSec(bed,0) << " ] seconds" << endl);
    }
    // Case 3: a dynamic single bed
    else if (m_nbTimeFrames>1 && m_nbBeds==1)
    {
      Cout("oImageDimensionsAndQuantification::SetAcquisitionTime() -> Dynamic single bed acquisition with following " << m_nbTimeFrames << " frame durations:" << endl);
      for (int fr=0; fr<m_nbTimeFrames; fr++)
        Cout("  --> Frame " << fr+1 << " with duration [ " << GetFrameTimeStartInSec(0,fr) << " : " << GetFrameTimeStopInSec(0,fr) << " ] seconds" << endl);
    }
    // Case 4: a dynamic multi beds
    else
    {
      Cout("oImageDimensionsAndQuantification::SetAcquisitionTime() -> Dynamic " << m_nbBeds << " beds acquistion with following " << m_nbTimeFrames << " frame durations:" << endl);
      for (int bed=0; bed<m_nbBeds; bed++)
      {
        Cout("  --> Bed " << bed+1 << " as following framing:" << endl);
        for (int fr=0; fr<m_nbTimeFrames; fr++)
          Cout("  |   Frame " << fr+1 << " with duration [ " << GetFrameTimeStartInSec(bed,fr) << " : " << GetFrameTimeStopInSec(bed,fr) << " ] seconds" << endl);
      }
    }
    // Correct for frame duration or not
    if (m_ignoreFdurCorrectionFlag) Cout("  --> Ignore frame duration correction" << endl);
    else Cout("  --> Correct for frame duration" << endl);
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::SetCalibrationFactor(int a_bed, FLTNB a_calibrationFactor)
{
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageDimensionsAndQuantification::SetCalibrationFactor() -> Object not initialized !" << endl);
    return 1;
  }
  // Ignore calibration factor if asked for
  if (m_ignoreCaliCorrectionFlag)
  {
    // Verbose
    if (m_verbose>=2 && a_bed==m_nbBeds-1) Cout("oImageDimensionsAndQuantification::SetCalibrationFactor() -> Ignore calibration factor correction" << endl);
    // Exit the function
    return 0;
  }
  // Check if calibration factor is strictly positive
  if (a_calibrationFactor<=0.)
  {
    Cerr("***** oImageDimensionsAndQuantification::SetCalibrationFactor() -> Provided calibration factor (" << a_calibrationFactor << ") is negative or null !" << endl);
    return 1;
  }
  // Affect quantification factor for all frames for this bed (even though the calibration factor should be the same
  // for all beds, we do not check it, it should be self consistent in the input files)
  if (!m_ignoreCaliCorrectionFlag)
  {
    for (int fr=0; fr<m_nbTimeFrames; fr++) for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
      m3p_quantificationFactors[a_bed][fr][g] *= a_calibrationFactor;
  }
  // Verbose
  if (m_verbose>=2 && a_bed==m_nbBeds-1)
    Cout("oImageDimensionsAndQuantification::SetCalibrationFactor() -> Correct for following calibration factor: " << a_calibrationFactor << endl);
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::SetDynamicSpecificQuantificationFactors(const string& a_quantificationFile)
{
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageDimensionsAndQuantification::SetDynamicSpecificQuantificationFactors() -> Object not initialized !" << endl);
    return 1;
  }

  // Provide the DynamicDataManager with a 2D array containing quantitative factors to update during initialization
  mp_DynamicDataManager->SetDynamicSpecificQuantificationFactors(m3p_quantificationFactors[0]);

  // Exit the function if no file provided
  if (a_quantificationFile=="") return 0;

  // Currently, it is not possible to reconstruction a multi-bed gated acquisition, so we throw an error in this function
  if (m_nbBeds>1)
  {
    Cerr("***** oImageDimensionsAndQuantification::SetDynamicSpecificQuantificationFactors() -> Multi-bed gated acquisitions cannot be reconstructed !" << endl);
    return 1;
  }

  // Verbose
  if (m_verbose>=2) Cout("oImageDimensionsAndQuantification::SetDynamicSpecificQuantificationFactors()-> Processing quantification file '" << a_quantificationFile << "'" << endl);
  
  // GET QUANTIFICATION CORRECTION FACTORS WITH RESPECT TO FRAMES or GATES
  
  // READ USER-DEFINED QUANTIFICATION FACTORS
  
  // TODO : For now we assume the configuration file holds :
  // TODO : - frame quantification factors : specific to each bed & frame
  // TODO : - gate quantification factors : specific to each frame and gates
  // TODO : Might require a new ReadDataASCIIFile function to recover this
  //
  // TODO : Additionally, we should be careful in documentation and explained which gate quantitfication factors are natively taken into account or not 
  
  
  // Quantification factors for each bed/frame/gate, intialized to 1.
  FLTNB*** pp_dynamic_quantification_factors = new FLTNB**[m_nbBeds];
  for (int bed=0 ; bed<m_nbBeds ; bed++)
  {
    pp_dynamic_quantification_factors[bed] = new FLTNB*[m_nbTimeFrames];
    for (int fr=0; fr<m_nbTimeFrames; fr++)
    {
      pp_dynamic_quantification_factors[bed][fr] = new FLTNB[m_nbRespGates*m_nbCardGates];
      for (int g=0; g<m_nbRespGates*m_nbCardGates; g++) pp_dynamic_quantification_factors[bed][fr][g] = 1.;
    }
  }

  // Build bed-related key words
  string *bed_name = new string[m_nbBeds + 1];
  
  for(int bed=0 ; bed<m_nbBeds ; bed++)
  {
    ostringstream oss( ostringstream::out );
    oss << "bed" << bed+1;
    bed_name[bed] = oss.str();
  }
    
  bed_name[m_nbBeds] = "eof";
 
  for (int bed=0 ; bed<m_nbBeds ; bed++)
  {
    int return_value = ReadDataASCIIFile(a_quantificationFile, "QUANTIFICATION_FACTORS", pp_dynamic_quantification_factors[bed], m_nbRespGates*m_nbCardGates, m_nbTimeFrames, KEYWORD_MANDATORY, bed_name[bed], bed_name[bed+1]);
    if (return_value<0) // string not found error
    {
      Cerr("***** oImageDimensionsAndQuantification::SetDynamicSpecificQuantificationFactors() -> Didn't found quantitative factors in file " << a_quantificationFile << " !" << endl);
      return 1;
    }
    else if(return_value == 1) // reading error
    {
      Cerr("***** oImageDimensionsAndQuantification::SetDynamicSpecificQuantificationFactors() -> An error occurred while trying to recover specific quantitative factors for frame !" << endl);
      return 1;
    }
    else if(return_value == 0) // correct reading
    {
      for(int bed=0 ; bed<m_nbBeds ; bed++)
        for (int fr=0; fr<m_nbTimeFrames; fr++)
          for (int g=0; g<m_nbRespGates*m_nbCardGates; g++)
            if (pp_dynamic_quantification_factors[bed][fr][g] <= 0)
            {
              Cerr("***** oImageDimensionsAndQuantification::SetDynamicSpecificQuantificationFactors() -> Provided quantification factor (" << pp_dynamic_quantification_factors[bed][fr][g] << ") is negative or null !" << endl);
              return 1;
            }
            else
            {
              m3p_quantificationFactors[bed][fr][g] *= pp_dynamic_quantification_factors[bed][fr][g];
            }
    }
  }

  // Delete temporary tabs
  for (int bed=0; bed<m_nbBeds; bed++)
  {
    for (int fr=0; fr<m_nbTimeFrames; fr++) delete pp_dynamic_quantification_factors[bed][fr];
    delete[] pp_dynamic_quantification_factors[bed];
  }
  delete[] pp_dynamic_quantification_factors;
  delete[] bed_name;

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB oImageDimensionsAndQuantification::GetQuantificationFactor(int a_bed, int a_frame, int a_respGate, int a_cardGate) 
{
  // Make sure that if respiratory or cardiac motion is enabled, we will recover the factor corresponding to the first image
  // TODO TM : perhaps more logical to initialize the oImageDimensions with the actual number of gates even when motion is enabled 
  // TODO TM : (i.e : m_nbRespGates or m_nbCardGates will be equal to 1, but the sensitivity image (list-mode) will require quantification factor for each gate
  // TODO TM : and add GetNbRespImages() functions anywhere in the code we need to know the number of resp/card images to be reconstructed
  // SS: If the motion correction is enabled, then all gates are pulled together in the optimization process, so there is no need for any quantification factor specific to each gate...
  if (m_nbRespGates == 1) a_respGate = 0;
  if (m_nbCardGates == 1) a_cardGate = 0;
  return m3p_quantificationFactors[a_bed][a_frame][a_respGate*m_nbCardGates+a_cardGate];
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::SetSPECTIsotope(int a_bed, const string& a_isotope)
{
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageDimensionsAndQuantification::SetSPECTIsotope() -> Object not initialized !" << endl);
    return 1;
  }

  // Not yet implemented

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::SetPETIsotope(int a_bed, const string& a_isotope)
{
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageDimensionsAndQuantification::SetPETIsotope() -> Object not initialized !" << endl);
    return 1;
  }

  // ------------------------------------------------------------------
  // Preliminary step is checking if we ignore isotope corrections
  // ------------------------------------------------------------------

  if (m_ignoreBratCorrectionFlag && m_ignoreDecaCorrectionFlag)
  {
    // Verbose
    if (m_verbose>=2 && a_bed==m_nbBeds-1)
      Cout("oImageDimensionsAndQuantification::SetPETIsotope() -> Ignore isotope dependent corrections" << endl);
    // Exit the function
    return 0;
  }

  // ------------------------------------------------------------------
  // Preliminary step is checking for unknown keyword
  // ------------------------------------------------------------------

  // Check if the isotope is named as unknown
  if (a_isotope=="UNKNOWN" || a_isotope=="Unknown" || a_isotope=="unknown")
  {
    // In this case we simply assume perfect branching ratio and infinite half-life
    if (m_verbose>=2 && a_bed==m_nbBeds-1)
      Cout("oImageDimensionsAndQuantification::SetPETIsotope() -> Un-specified isotope; no decay nor branching ratio correction" << endl);
    // Exit the function
    return 0;
  }

  // ------------------------------------------------------------------
  // First step is open the isotopes file and find the provided isotope
  // ------------------------------------------------------------------

  // Get oOutputManager instance and config directory
  sOutputManager* p_output_manager = sOutputManager::GetInstance();
  string config_dir = p_output_manager->GetPathToConfigDir();

  // Open the isotope file based on config directory
  string file_name = config_dir + "/misc/isotopes_pet.txt";
  ifstream fin(file_name.c_str());

  if (!fin)
  {
    Cerr("***** oImageDimensionsAndQuantification::SetPETIsotope() -> Failed to open PET isotopes data file '" << file_name << "' !" << endl);
    return 1;
  }

  // Loop on lines to find the isotope
  int line_max_size = 10240;
  char *line = new char[line_max_size];
  FLTNB half_life = -1.;
  FLTNB branching_ratio = -1.;
  bool found_it = false;
  fin.getline(line,line_max_size);
  while (!fin.eof())
  {
    // For the word position
    size_t found_position;
    // Cast the line into a string for convenience
    string test = (string)line;
    // Jump to next line if we find a # character as the first character
    if ((found_position=test.find("#"))==0)
    {
      // Read a new line before continuing
      fin.getline(line,line_max_size);
      continue;
    }
    // Check if we see the isotope name
    found_position = test.find(a_isotope);
    if (found_position!=string::npos)
    {
      // Each line is organised as follows: (some spaces/tabs) isotope_name (some spaces/tabs) half_life (some spaces/tabs) branching_ratio
      // We first remove the isotope name from the line
      test = test.substr(found_position+1+a_isotope.length());
      // Then we convert it into a string stream to ease the reading of both numbers
      istringstream fstr(test);
      fstr >> half_life >> branching_ratio;
      // We found it
      found_it = true;
      break;
    }
    // Read a new line
    fin.getline(line,line_max_size);
  }

  delete[] line;
  // Close file
  fin.close();

  // Check if we found it or not
  if (found_it)
  {
    // Check rationality of values
    if (branching_ratio<=0. || branching_ratio>1.)
    {
      Cerr("***** oImageDimensionsAndQuantification::SetPETIsotope() -> Branching ratio (" << branching_ratio << ") is not in the ]0:1] range !" << endl);
      return 1;
    }
    // Verbose
    if (m_verbose>=2 && a_bed==m_nbBeds-1)
    {
      // If negative half-life, then we consider it is infinite (no half-life)
      if (half_life<=0.)
        Cout("oImageDimensionsAndQuantification::SetPETIsotope() -> Isotope " << a_isotope << " has infinite half life and " << branching_ratio << " branching ratio" << endl);
      else
        Cout("oImageDimensionsAndQuantification::SetPETIsotope() -> Isotope " << a_isotope << " has " << half_life << " seconds half life and " << branching_ratio << " branching ratio" << endl);
    }
  }
  else
  {
    // Throw error message
    Cerr("***** oImageDimensionsAndQuantification::SetPETIsotope() -> Did not find " << a_isotope << " isotope in the PET isotope data file, please add it !" << endl);
    return 1;
  }

  // ------------------------------------------------------------------------
  // Second step is applying branching ratio and decay wrt to frame durations
  // ------------------------------------------------------------------------

  // Branching ratio
  if (!m_ignoreBratCorrectionFlag)
  {
    // Apply correction
    for (int fr=0; fr<m_nbTimeFrames; fr++)
    {
      for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
        m3p_quantificationFactors[a_bed][fr][g] /= branching_ratio;
    }
    // Verbose
    if (m_verbose>=2 && a_bed==m_nbBeds-1) Cout("  --> Correct for branching ratio" << endl);
  }
  // Verbose
  else if (m_verbose>=2 && a_bed==m_nbBeds-1) Cout("  --> Ignore branching ratio correction" << endl);

  // Apply decay factors if half life is not negative (meaning no half-life)
  if (half_life>0.)
  {
    // We correct for it
    if (!m_ignoreDecaCorrectionFlag)
    {
      // Apply correction
      for (int fr=0; fr<m_nbTimeFrames; fr++)
      {
        m_lambda = log(2.0)/half_life;
        for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
        {
          long double dstart = m_lambda*GetFrameTimeStartInSec(a_bed,fr);
          long double dacq = m_lambda*GetFrameDurationInSec(a_bed,fr);
          // Time start decay correction
          m3p_quantificationFactors[a_bed][fr][g] *= exp(dstart);
          // Frame duration decay correction
          m3p_quantificationFactors[a_bed][fr][g] *= dacq/(1.0-exp(-dacq));
          
          /* Time start decay correction
          m3p_quantificationFactors[a_bed][fr][g] *= exp(lambda*GetFrameTimeStartInSec(a_bed,fr));
          // Frame duration decay correction
          m3p_quantificationFactors[a_bed][fr][g] *= lambda*GetFrameDurationInSec(a_bed,fr)/(1.0-exp(-lambda*(GetFrameDurationInSec(a_bed,fr))));*/
        }
      }
      // Verbose
      if (m_verbose>=2 && a_bed==m_nbBeds-1) Cout("  --> Correct for half-life" << endl);
    }
    // Verbose
    else if (m_verbose>=2 && a_bed==m_nbBeds-1) Cout("  --> Ignore half-life correction" << endl);
  }
  
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::InitDynamicData( string a_pathTo4DDataSplittingFile,
                                                        int a_respMotionCorrectionFlag,
                                                        int a_cardMotionCorrectionFlag,
                                                        int a_invMotionCorrectionFlag,
                                                        int a_nbRespGates, int a_nbCardGates )
{
  if (m_verbose>=5) Cout("oImageDimensionsAndQuantification::InitDynamicData()" << endl);
  mp_DynamicDataManager->SetVerbose(m_verbose);
  mp_DynamicDataManager->SetImageDimensionsAndQuantification(this);
  if (mp_DynamicDataManager->InitDynamicData( a_nbRespGates, a_nbCardGates, a_pathTo4DDataSplittingFile,
                                              a_respMotionCorrectionFlag, a_cardMotionCorrectionFlag, a_invMotionCorrectionFlag ))
  {
    Cerr("***** oImageDimensionsAndQuantification::InitDynamicData() -> A problem occurred while initializing the dynamic data from dynamic data manager !" << endl);
    return 1;
  }
  
  // Set type of dynamic reconstruction
  if (a_respMotionCorrectionFlag || a_cardMotionCorrectionFlag)
    m_dynRecoTypeFlag = DYN_RECO_MCGATING;
  else if (a_nbRespGates>1 || a_nbCardGates>1)
    m_dynRecoTypeFlag = DYN_RECO_GATING;
  else if (a_invMotionCorrectionFlag)
    m_dynRecoTypeFlag = DYN_RECO_IPMC;
  else if (m_nbTimeFrames>1)
    m_dynRecoTypeFlag = DYN_RECO_FRAMING;
    
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageDimensionsAndQuantification::CheckDynamicParameters(int64_t a_nbEvents)
{
  if (m_verbose>=5) Cout("oImageDimensionsAndQuantification::CheckDynamicParameters()" << endl);
  return mp_DynamicDataManager->CheckParameters(a_nbEvents);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
