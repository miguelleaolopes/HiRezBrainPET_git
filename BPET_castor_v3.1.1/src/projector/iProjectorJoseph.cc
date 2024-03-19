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
  \ingroup  projector
  \brief    Implementation of class iProjectorJoseph
*/

#include "iProjectorJoseph.hh"
#include "sOutputManager.hh"

#include <cmath>
#ifdef _WIN32
// Avoid compilation errors due to mix up between std::min()/max() and
// min max macros
#undef min
#undef max
#endif

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorJoseph::iProjectorJoseph() : vProjector()
{
  // This projector is not compatible with SPECT attenuation correction because
  // the voxels contributing to the line are not strictly ordered with respect to
  // their distance to point2 (due to interpolations at each plane crossed)
  m_compatibleWithSPECTAttenuationCorrection = false;
  // This projector is compatible with compression as it works only with the
  // cartesian coordinates of 2 end points of a line
  m_compatibleWithCompression = true;
  // Default pointers and parameters
  mp_boundariesX = nullptr;
  mp_boundariesY = nullptr;
  mp_boundariesZ = nullptr;
  mp_maskPad = nullptr;
  m_toleranceX = 0.;
  m_toleranceY = 0.;
  m_toleranceZ = 0.;
  m_tolerance_fctr = 1.e-6;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorJoseph::~iProjectorJoseph()
{
  if ( mp_boundariesX )
  {
    delete[] mp_boundariesX;
    mp_boundariesX = nullptr;
  }

  if ( mp_boundariesY )
  {
    delete[] mp_boundariesY;
    mp_boundariesY = nullptr;
  }

  if ( mp_boundariesZ )
  {
    delete[] mp_boundariesZ;
    mp_boundariesZ = nullptr;
  }

  if ( mp_maskPad )
  {
    delete[] mp_maskPad;
    mp_maskPad = nullptr;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorJoseph::ReadConfigurationFile(const string& a_configurationFile)
{
  // No options for joseph
  ;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorJoseph::ReadOptionsList(const string& a_optionsList)
{
  // Currently only one option: Tolerance with respect to voxel sizes in each dimensions
  uint16_t option[1];
  
  // Check tolerance value
  if (ReadStringOption(a_optionsList, option, 1, ",", "Tolerance factor for plane index computation "))
  {
    Cerr("***** iProjectorJoseph::ReadConfigurationFile() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  
  // Check tolerance factor
  if (option[0] < 2.)
  {
    Cerr("***** iProjectorJoseph::ReadOptionsList() -> Tolerance factor parameter must be  >2 (tolerance factor < 10^-2)" << endl);
    Cerr("                                          -> Provided parameter: "<< option[0] << endl);
    return 1;
  }

  m_tolerance_fctr = pow(10 , -(HPFLTNB)option[0]);

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iProjectorJoseph::ShowHelpSpecific()
{
  cout << "This projector is a line projector that uses linear interpolation between pixels." << endl;
  cout << "It is implemented from the following published paper:" << endl;
  cout << "P. M. Joseph, \"An improved algorithm for reprojecting rays through pixel images\", IEEE Trans. Med. Imaging, vol. 1, pp. 192-6, 1982." << endl;
  cout << "There is 1 optional parameter for this projector:" << endl;
  cout << "  tolerance_factor: x. -> 'x' must be an unsigned integer, and represent the tolerance precision factor for plane index computation (10^(-x))." << endl;
  cout << "                       -> Default: 'x' = 6" << endl;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorJoseph::CheckSpecificParameters()
{
  if ( m_tolerance_fctr >  0.01 
    || m_tolerance_fctr <= 0 )
  {
    Cerr("***** iProjectorJoseph::CheckSpecificParameters() -> Tolerance factor must be defined in the specific interval : ]0 ; 1.e-2]" << endl);
    Cerr("*****                                             -> Default value is 1.e-6" << endl);
    Cerr("*****                                             -> Current value is " << m_tolerance_fctr << endl);
    return 1;
  }
  
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorJoseph::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=2) Cout("iProjectorJoseph::InitializeSpecific() -> Use Joseph projector" << endl);

  // Allocate and compute boundaries (grid through voxel centers)
  mp_boundariesX = new HPFLTNB[ mp_nbVox[ 0 ] + 2 ];
  for( INTNB i = 0; i < mp_nbVox[ 0 ] + 2; ++i ) mp_boundariesX[ i ] = -((HPFLTNB)(mp_halfFOV[ 0 ])) + ((HPFLTNB)(mp_sizeVox[ 0 ])) * ( ((HPFLTNB)i) - 0.5 );
  mp_boundariesY = new HPFLTNB[ mp_nbVox[ 1 ] + 2 ];
  for( INTNB i = 0; i < mp_nbVox[ 1 ] + 2; ++i ) mp_boundariesY[ i ] = -((HPFLTNB)(mp_halfFOV[ 1 ])) + ((HPFLTNB)(mp_sizeVox[ 1 ])) * ( ((HPFLTNB)i) - 0.5 );
  mp_boundariesZ = new HPFLTNB[ mp_nbVox[ 2 ] + 2 ];
  for( INTNB i = 0; i < mp_nbVox[ 2 ] + 2; ++i ) mp_boundariesZ[ i ] = -((HPFLTNB)(mp_halfFOV[ 2 ])) + ((HPFLTNB)(mp_sizeVox[ 2 ])) * ( ((HPFLTNB)i) - 0.5 );

  // Allocating the mask buffer for the padded image space
  INTNB nElts = ( mp_nbVox[ 0 ] + 2 ) * ( mp_nbVox[ 1 ] + 2 ) * ( mp_nbVox[ 2 ] + 2 );
  mp_maskPad = new uint8_t[ nElts ];
  ::memset( mp_maskPad, 0, sizeof( uint8_t ) * nElts );
  for( INTNB k = 1; k < mp_nbVox[ 2 ] + 1; ++k )
  {
    for( INTNB j = 1; j < mp_nbVox[ 1 ] + 1; ++j )
    {
      for( INTNB i = 1; i < mp_nbVox[ 0 ] + 1; ++i )
      {
        mp_maskPad[ i + j * ( ( mp_nbVox[ 0 ] + 2 ) ) + k * ( mp_nbVox[ 0 ] + 2 ) * ( mp_nbVox[ 1 ] + 2 ) ] = 1;
      }
    }
  }
  
  // Set the tolerance with respect to voxel sizes in each dimensions
  m_toleranceX = ((HPFLTNB)(mp_sizeVox[0])) * m_tolerance_fctr;
  m_toleranceY = ((HPFLTNB)(mp_sizeVox[1])) * m_tolerance_fctr;
  m_toleranceZ = ((HPFLTNB)(mp_sizeVox[2])) * m_tolerance_fctr;

  // Setting the bounds
  m_boundX = (-mp_halfFOV[ 0 ]) - mp_sizeVox[ 0 ] * 0.5;
  m_boundY = (-mp_halfFOV[ 1 ]) - mp_sizeVox[ 1 ] * 0.5;
  m_boundZ = (-mp_halfFOV[ 2 ]) - mp_sizeVox[ 2 ] * 0.5;

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

INTNB iProjectorJoseph::EstimateMaxNumberOfVoxelsPerLine()
{
  // Find the maximum number of voxels along a given dimension
  INTNB max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxX();
  if (mp_ImageDimensionsAndQuantification->GetNbVoxY()>max_nb_voxels_in_dimension) max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxY();
  if (mp_ImageDimensionsAndQuantification->GetNbVoxZ()>max_nb_voxels_in_dimension) max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxZ();
  // We should have at most 4 voxels in a given plane, so multiply by 4
  max_nb_voxels_in_dimension *= 4;
  // Return the value
  return max_nb_voxels_in_dimension;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorJoseph::ProjectWithoutTOF(int a_direction, oProjectionLine* ap_ProjectionLine )
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorJoseph::ProjectWithoutTOF() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorJoseph::Project without TOF -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // Get event positions
  FLTNB* event1 = ap_ProjectionLine->GetPosition1();
  FLTNB* event2 = ap_ProjectionLine->GetPosition2();

  // Distance between point event1 and event2
  HPFLTNB const r[ 3 ] = {
    event2[ 0 ] - event1[ 0 ],
    event2[ 1 ] - event1[ 1 ],
    event2[ 2 ] - event1[ 2 ]
  };

  // Square of r
  HPFLTNB const r2[ 3 ] = {
    r[ 0 ] * r[ 0 ],
    r[ 1 ] * r[ 1 ],
    r[ 2 ] * r[ 2 ]
  };

  // Find the first and last intersecting plane in X using the parametric
  // values alphaMin and alphaMax
  HPFLTNB alphaMin = 0.0, alphaMax = 1.0;

  // Variables for Joseph
  HPFLTNB pos[ 3 ] = { 0.0, 0.0, 0.0 }; // Position of point in image space
  HPFLTNB wfl[ 2 ] = { 0.0, 0.0 }; // Weight floor
  HPFLTNB wcl[ 2 ] = { 0.0, 0.0 }; // Weight ceil
  HPFLTNB w[ 4 ] = { 0.0, 0.0, 0.0, 0.0 }; // Interpolation weight
  int16_t index[ 2 ] = { 0, 0 }; // index in the image space
  int32_t finalIdx = 0; // final index
  int8_t limitX1 = 1; int8_t limitX2 = 1;
  int8_t limitY1 = 1; int8_t limitY2 = 1;
  int8_t limitZ1 = 1; int8_t limitZ2 = 1;

  // Condition on the largest component of r
  if( ::fabs( r[ 0 ] ) > ::fabs( r[ 1 ] ) )
  {
    // Computing the parametric values
    // For the X-axis
    // We stay in the limit of the image space
    HPFLTNB const alphaX_0 = ( -mp_halfFOV[ 0 ] - event1[ 0 ] ) / r[ 0 ];
    HPFLTNB const alphaX_1 = (  mp_halfFOV[ 0 ] - event1[ 0 ] ) / r[ 0 ];
    alphaMin = std::max( alphaMin, std::min( alphaX_0, alphaX_1 ) );
    alphaMax = std::min( alphaMax, std::max( alphaX_0, alphaX_1 ) );

    // For the Y-axis
    // We introduce 1 voxel size around Y-axis
    // Make a shift on first and last plane (like an offset)
    if( r[ 1 ] != 0.0 )
    {
      HPFLTNB const alphaY_0 = ( -mp_halfFOV[ 1 ] - mp_sizeVox[ 1 ] * 0.5 - event1[ 1 ] )
        / r[ 1 ];
      HPFLTNB const alphaY_1 = ( mp_halfFOV[ 1 ] + mp_sizeVox[ 1 ] * 0.5 - event1[ 1 ] )
        / r[ 1 ];
      alphaMin = std::max( alphaMin, std::min( alphaY_0, alphaY_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaY_0, alphaY_1 ) );
    }

    // For the Z-axis 
    // We introduce 1 voxel size around Z-axis
    // Make a shift on first and last plane (like an offset)
    if( r[ 2 ] != 0.0 )
    {
      HPFLTNB const alphaZ_0 = ( -mp_halfFOV[ 2 ] - mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      HPFLTNB const alphaZ_1 = ( mp_halfFOV[ 2 ] + mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      alphaMin = std::max( alphaMin, std::min( alphaZ_0, alphaZ_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaZ_0, alphaZ_1 ) );
    }

    if( alphaMax <= alphaMin ) return 0;

    if( r[ 1 ] == 0.0 ) if( event1[ 1 ] < -mp_halfFOV[ 1 ] || event1[ 1 ] > mp_halfFOV[ 1 ] ) return 0;
    if( r[ 2 ] == 0.0 ) if( event1[ 2 ] < -mp_halfFOV[ 2 ] || event1[ 2 ] > mp_halfFOV[ 2 ] ) return 0;

    // Computing weight of normalization
    HPFLTNB const factor( ::fabs( r[ 0 ] )
      / ::sqrt( r2[ 0 ] + r2[ 1 ] + r2[ 2 ] ) );
    HPFLTNB const weight( mp_sizeVox[ 0 ] / factor );

    // Computing the increment
    HPFLTNB const ri[ 3 ] = {
      1.0, // r[ 0 ] / r[ 0 ]
      r[ 1 ] / r[ 0 ],
      r[ 2 ] / r[ 0 ]
    };

    // Computing the first and the last plane
    int iMin = 0, iMax = 0;
    if( r[ 0 ] > 0.0 )
    {
      iMin = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alphaMin * r[ 0 ] - event1[ 0 ] ) / mp_sizeVox[ 0 ] );
      iMax = ::floor( m_toleranceX + ( event1[ 0 ] + alphaMax * r[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }
    else if( r[ 0 ] < 0.0 )
    {
      iMin = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alphaMax * r[ 0 ] - event1[ 0 ] ) / mp_sizeVox[ 0 ] );
      iMax = ::floor( m_toleranceX + ( event1[ 0 ] + alphaMin * r[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = -mp_halfFOV[ 0 ] + ( mp_sizeVox[ 0 ] * 0.5 ) - event1[ 0 ];

    // Loop over all the crossing planes
    for( int i = iMin; i < iMax; ++i )
    {
      // Computing position on crossed plane in term of grid spacing
      HPFLTNB const step = offset + i * mp_sizeVox[ 0 ];
      pos[ 1 ] = event1[ 1 ] + step * ri[ 1 ];
      pos[ 2 ] = event1[ 2 ] + step * ri[ 2 ];

      // Find the index in the image
      index[ 0 ] = ::floor( m_toleranceY + ( pos[ 1 ] - m_boundY ) / mp_sizeVox[ 1 ] );
      index[ 1 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - m_boundZ ) / mp_sizeVox[ 2 ] );

      // Storing values and indices
      finalIdx = i + ( index[ 0 ] - 1 ) * mp_nbVox[ 0 ] + ( index[ 1 ] - 1 ) * m_nbVoxXY;

      // Bilinear interpolation component (using floor)
      wfl[ 0 ] = pos[ 1 ] - ( m_boundY + index[ 0 ] * mp_sizeVox[ 1 ] );
      wfl[ 1 ] = pos[ 2 ] - ( m_boundZ + index[ 1 ] * mp_sizeVox[ 2 ] );
      wfl[ 0 ] /= mp_sizeVox[ 1 ];
      wfl[ 1 ] /= mp_sizeVox[ 2 ];

      // Bilinear interpolation component (using ceil)
      wcl[ 0 ] = 1.0 - wfl[ 0 ];
      wcl[ 1 ] = 1.0 - wfl[ 1 ];

      // Final coefficients
      w[ 0 ] = wcl[ 0 ] * wcl[ 1 ];
      w[ 1 ] = wcl[ 0 ] * wfl[ 1 ];
      w[ 2 ] = wfl[ 0 ] * wfl[ 1 ];
      w[ 3 ] = wfl[ 0 ] * wcl[ 1 ];

      // Check Y bounds
      if( index[ 0 ] <= 0 ) limitY1 = 0;
      if( index[ 0 ] >= mp_nbVox[ 1 ] ) limitY2 = 0;

      // Check Z bounds
      if( index[ 1 ] <= 0  ) limitZ1 = 0;
      if( index[ 1 ] >= mp_nbVox[ 2 ] ) limitZ2 = 0;

      if (!m_hasMask || (limitY1 && limitZ1 && index[ 0 ] <= mp_nbVox[ 1 ] && index[ 1 ] <= mp_nbVox[ 2 ] && mp_mask[finalIdx]))
      {
        ap_ProjectionLine->AddVoxel(a_direction, finalIdx * limitY1 * limitZ1, w[ 0 ] * weight * limitY1 * limitZ1);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + m_nbVoxXY ) * limitY1 * limitZ2, w[ 1 ] * weight * limitY1 * limitZ2);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + mp_nbVox[ 0 ] + m_nbVoxXY ) * limitY2 * limitZ2, w[ 2 ] * weight * limitY2 * limitZ2);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + mp_nbVox[ 0 ] ) * limitY2 * limitZ1, w[ 3 ] * weight * limitY2 * limitZ1);
      }
      limitY1 = 1; limitY2 = 1;
      limitZ1 = 1; limitZ2 = 1;
    }
  }
  else
  {
    // Computing the parametric values
    // For the X-axis
    // We introduce 1 voxel size around Y-axis
    if( r[ 0 ] != 0.0 )
    {
      HPFLTNB const alphaX_0 = ( -mp_halfFOV[ 0 ] - mp_sizeVox[ 0 ] * 0.5 - event1[ 0 ] )
        / r[ 0 ];
      HPFLTNB const alphaX_1 = ( mp_halfFOV[ 0 ] + mp_sizeVox[ 0 ] * 0.5 - event1[ 0 ] )
        / r[ 0 ];
      alphaMin = std::max( alphaMin, std::min( alphaX_0, alphaX_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaX_0, alphaX_1 ) );
    }

    // For the Y-axis
    // Make a shift on first and last plane (like an offset)
    // We stay in the limit of the image space
    HPFLTNB const alphaY_0 = ( -mp_halfFOV[ 1 ] - event1[ 1 ] ) / r[ 1 ];
    HPFLTNB const alphaY_1 = ( mp_halfFOV[ 1 ] - event1[ 1 ] ) / r[ 1 ];
    alphaMin = std::max( alphaMin, std::min( alphaY_0, alphaY_1 ) );
    alphaMax = std::min( alphaMax, std::max( alphaY_0, alphaY_1 ) );

    // For the Z-axis 
    // We introduce 1 voxel size around Z-axis
    // Make a shift on first and last plane (like an offset)
    if( r[ 2 ] != 0.0 )
    {
      HPFLTNB const alphaZ_0 = ( -mp_halfFOV[ 2 ] - mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      HPFLTNB const alphaZ_1 = ( mp_halfFOV[ 2 ] + mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      alphaMin = std::max( alphaMin, std::min( alphaZ_0, alphaZ_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaZ_0, alphaZ_1 ) );
    }

    if( alphaMax <= alphaMin ) return 0;

    if( r[ 0 ] == 0.0 ) if( event1[ 0 ] < -mp_halfFOV[ 0 ] || event1[ 0 ] > mp_halfFOV[ 0 ] ) return 0;
    if( r[ 2 ] == 0.0 ) if( event1[ 2 ] < -mp_halfFOV[ 2 ] || event1[ 2 ] > mp_halfFOV[ 2 ] ) return 0;

    // Computing weight of normalization
    HPFLTNB const factor( ::fabs( r[ 1 ] ) / ::sqrt( r2[ 0 ] + r2[ 1 ] + r2[ 2 ] ) );
    HPFLTNB const weight( mp_sizeVox[ 1 ] / factor );

    // Computing the increment
    HPFLTNB const ri[ 3 ] = {
      r[ 0 ] / r[ 1 ],
      1.0, // r[ 1 ] / r[ 1 ]
      r[ 2 ] / r[ 1 ]
    };

    // Computing the first and the last plane
    int jMin = 0, jMax = 0;
    if( r[ 1 ] > 0.0 )
    {
      jMin = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alphaMin * r[ 1 ] - event1[ 1 ] ) / mp_sizeVox[ 1 ] );
      jMax = ::floor( m_toleranceY + ( event1[ 1 ] + alphaMax * r[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }
    else if( r[ 1 ] < 0.0 )
    {
      jMin = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alphaMax * r[ 1 ] - event1[ 1 ] ) / mp_sizeVox[ 1 ] );
      jMax = ::floor( m_toleranceY + ( event1[ 1 ] + alphaMin * r[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = -mp_halfFOV[ 1 ] + ( mp_sizeVox[ 1 ] * 0.5 ) - event1[ 1 ];

    // Loop over all the crossing planes
    for( int j = jMin; j < jMax; ++j )
    {
      // Computing position on crossed plane in term of grid spacing
      HPFLTNB const step = offset + j * mp_sizeVox[ 1 ];
      pos[ 0 ] = event1[ 0 ] + step * ri[ 0 ];
      pos[ 2 ] = event1[ 2 ] + step * ri[ 2 ];

      // Find the index in the image
      index[ 0 ] = ::floor( m_toleranceX + ( pos[ 0 ] - m_boundX ) / mp_sizeVox[ 0 ] );
      index[ 1 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - m_boundZ ) / mp_sizeVox[ 2 ] );

      // Storing values and indices
      finalIdx = ( index[ 0 ] - 1 ) + j * mp_nbVox[ 0 ] + ( index[ 1 ] - 1 ) * m_nbVoxXY;

      // Bilinear interpolation component (using floor)
      wfl[ 0 ] = pos[ 0 ] - ( m_boundX + index[ 0 ] * mp_sizeVox[ 0 ] );
      wfl[ 1 ] = pos[ 2 ] - ( m_boundZ + index[ 1 ] * mp_sizeVox[ 2 ] );
      wfl[ 0 ] /= mp_sizeVox[ 0 ];
      wfl[ 1 ] /= mp_sizeVox[ 2 ];

      // Bilinear interpolation component (using ceil)
      wcl[ 0 ] = 1.0 - wfl[ 0 ];
      wcl[ 1 ] = 1.0 - wfl[ 1 ];

      // Final coefficients
      w[ 0 ] = wcl[ 0 ] * wcl[ 1 ];
      w[ 1 ] = wcl[ 0 ] * wfl[ 1 ];
      w[ 2 ] = wfl[ 0 ] * wfl[ 1 ];
      w[ 3 ] = wfl[ 0 ] * wcl[ 1 ];

      // Check X bounds
      if( index[ 0 ] <= 0 ) limitX1 = 0;
      if( index[ 0 ] >= mp_nbVox[ 0 ] ) limitX2 = 0;

      // Check Z bounds
      if( index[ 1 ] <= 0 ) limitZ1 = 0;
      if( index[ 1 ] >= mp_nbVox[ 2 ] ) limitZ2 = 0;

      if (!m_hasMask || (limitX1 && limitZ1 && index[ 0 ] <= mp_nbVox[ 0 ] && index[ 1 ] <= mp_nbVox[ 2 ] && mp_mask[finalIdx]))
      {
        ap_ProjectionLine->AddVoxel(a_direction, finalIdx * limitX1 * limitZ1, w[ 0 ] * weight * limitX1 * limitZ1);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + m_nbVoxXY ) * limitX1 * limitZ2, w[ 1 ] * weight * limitX1 * limitZ2);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + 1 + m_nbVoxXY ) * limitX2 * limitZ2, w[ 2 ] * weight * limitX2 * limitZ2);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + 1 ) * limitX2 * limitZ1, w[ 3 ] * weight * limitX2 * limitZ1);
      }

      limitX1 = 1; limitX2 = 1;
      limitZ1 = 1; limitZ2 = 1;
    }
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorJoseph::ProjectTOFListmode(int a_direction, oProjectionLine* ap_ProjectionLine)
{

  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorJoseph::ProjectTOFListmode() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorJoseph::Project with TOF measurement -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // Get event positions
  FLTNB* event1 = ap_ProjectionLine->GetPosition1();
  FLTNB* event2 = ap_ProjectionLine->GetPosition2();

  // Distance between point event1 and event2
  HPFLTNB const r[ 3 ] = {
    event2[ 0 ] - event1[ 0 ],
    event2[ 1 ] - event1[ 1 ],
    event2[ 2 ] - event1[ 2 ]
  };

/*
  // Square of r
  HPFLTNB const r2[ 3 ] = {
    r[ 0 ] * r[ 0 ],
    r[ 1 ] * r[ 1 ],
    r[ 2 ] * r[ 2 ]
  };
*/

  // LOR length
  HPFLTNB lor_length = ap_ProjectionLine->GetLength();

  // Get TOF info
  FLTNB tof_measurement = ap_ProjectionLine->GetTOFMeasurementInPs();

  // convert delta time into delta length
  HPFLTNB tof_delta = tof_measurement * SPEED_OF_LIGHT_IN_MM_PER_PS * 0.5;

  // TOF standard deviation and truncation
  HPFLTNB tof_sigma = m_TOFResolutionInMm / TWO_SQRT_TWO_LN_2;
  HPFLTNB tof_sigma_sqrt2 = sqrt(2.)*tof_sigma;
  HPFLTNB tof_half_span = 0.;

  if (m_TOFWeightingFcnPrecomputedFlag) tof_half_span = ((HPFLTNB)m_TOFWeightingFcnNbSamples)/(2.*m_TOFPrecomputedSamplingFactor);
  else if (m_TOFBinProperProcessingFlag) tof_half_span = tof_sigma * m_TOFNbSigmas + 0.5 * m_TOFBinSizeInMm;
  else tof_half_span = tof_sigma * m_TOFNbSigmas;

  // distance between the first event1 and the TOF measurement
  HPFLTNB lor_tof_center = lor_length * 0.5 + tof_delta;

  // coordinates of the lower edge of the first voxel falling inside the truncated TOF distribution centered at the TOF measurement
  HPFLTNB tof_edge_low[] = {0,0,0};
  // coordinate of the upper edge of the last voxel falling inside the truncated TOF distribution centered at the TOF measurement
  HPFLTNB tof_edge_high[] = {0,0,0};
  HPFLTNB tof_center;
  INTNB tof_index;
  
  // low/high voxel edges (in absolute coordinates) for truncated TOF
  for (int ax=0;ax<3;ax++)
  {
    // absolute coordinate along each axis of the center of the TOF distribution
    tof_center = event1[ax] +  lor_tof_center * r[ax] / lor_length;

    // absolute coordinate along each axis of the lowest voxel edge spanned by the TOF distribution, limited by the lowest FOV edge
    tof_edge_low[ax] = tof_center - tof_half_span * fabs(r[ax]) / lor_length;
    tof_index = max( (INTNB)::floor( (tof_edge_low[ax] - (-mp_halfFOV[ax])) / mp_sizeVox[ax] ), (INTNB)0);
    // if low TOF edge above the highest FOV edge, return empty line
    if (tof_index>mp_nbVox[ax]-1) return 0;
    tof_edge_low[ax] = (HPFLTNB)tof_index *  mp_sizeVox[ax] - mp_halfFOV[ax];

    // absolute coordinate along each axis of the highest voxel edge spanned by the TOF distribution, limited by the highest FOV edge
    tof_edge_high[ax] = tof_center + tof_half_span * fabs(r[ax]) / lor_length;
    tof_index = min( (INTNB)::floor( (tof_edge_high[ax] - (-mp_halfFOV[ax])) / mp_sizeVox[ax] ), mp_nbVox[ax]-1);
    // if high TOF edge below the lowest FOV edge, return empty line
    if (tof_index<0) return 0;
    tof_edge_high[ax] = (HPFLTNB)(tof_index+1) * mp_sizeVox[ax] - mp_halfFOV[ax];    
  }

  // Find the first and last intersecting plane in X using the parametric
  // values alphaMin and alphaMax
  HPFLTNB alphaMin = 0.0, alphaMax = 1.0;

  // Variables for Joseph
  HPFLTNB pos[ 3 ] = { 0.0, 0.0, 0.0 }; // Position of point in image space
  HPFLTNB wfl[ 2 ] = { 0.0, 0.0 }; // Weight floor
  HPFLTNB wcl[ 2 ] = { 0.0, 0.0 }; // Weight ceil
  HPFLTNB w[ 4 ] = { 0.0, 0.0, 0.0, 0.0 }; // Interpolation weight
  int16_t index[ 2 ] = { 0, 0 }; // index in the image space
  int32_t finalIdx = 0; // final index
  int8_t limitX1 = 1; int8_t limitX2 = 1;
  int8_t limitY1 = 1; int8_t limitY2 = 1;
  int8_t limitZ1 = 1; int8_t limitZ2 = 1;

  // Condition on the largest component of r
  if( ::fabs( r[ 0 ] ) > ::fabs( r[ 1 ] ) )
  {
    // Computing the parametric values
    // For the X-axis
    // We stay in the limit of the image space
    HPFLTNB const alphaX_0 = ( tof_edge_low[ 0 ] - event1[ 0 ] ) / r[ 0 ];
    HPFLTNB const alphaX_1 = (  tof_edge_high[ 0 ] - event1[ 0 ] ) / r[ 0 ];
    alphaMin = std::max( alphaMin, std::min( alphaX_0, alphaX_1 ) );
    alphaMax = std::min( alphaMax, std::max( alphaX_0, alphaX_1 ) );

    // For the Y-axis
    // We introduce 1 voxel size around Y-axis
    // Make a shift on first and last plane (like an offset)
    if( r[ 1 ] != 0.0 )
    {
      HPFLTNB const alphaY_0 = ( tof_edge_low[ 1 ] - mp_sizeVox[ 1 ] * 0.5 - event1[ 1 ] )
        / r[ 1 ];
      HPFLTNB const alphaY_1 = ( tof_edge_high[ 1 ] + mp_sizeVox[ 1 ] * 0.5 - event1[ 1 ] )
        / r[ 1 ];
      alphaMin = std::max( alphaMin, std::min( alphaY_0, alphaY_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaY_0, alphaY_1 ) );
    }

    // For the Z-axis 
    // We introduce 1 voxel size around Z-axis
    // Make a shift on first and last plane (like an offset)
    if( r[ 2 ] != 0.0 )
    {
      HPFLTNB const alphaZ_0 = ( tof_edge_low[ 2 ] - mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      HPFLTNB const alphaZ_1 = ( tof_edge_high[ 2 ] + mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      alphaMin = std::max( alphaMin, std::min( alphaZ_0, alphaZ_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaZ_0, alphaZ_1 ) );
    }

    if( alphaMax <= alphaMin ) return 0;

    if( r[ 1 ] == 0.0 ) if( event1[ 1 ] < -mp_halfFOV[ 1 ] || event1[ 1 ] > mp_halfFOV[ 1 ] ) return 0;
    if( r[ 2 ] == 0.0 ) if( event1[ 2 ] < -mp_halfFOV[ 2 ] || event1[ 2 ] > mp_halfFOV[ 2 ] ) return 0;

    // Computing weight of normalization
    HPFLTNB const factor( ::fabs( r[ 0 ] ) / lor_length );
    HPFLTNB const factor_for_tof( lor_length / r[ 0 ]) ;
    HPFLTNB const weight( mp_sizeVox[ 0 ] / factor );

    // Computing the increment
    HPFLTNB const ri[ 3 ] = {
      1.0, // r[ 0 ] / r[ 0 ]
      r[ 1 ] / r[ 0 ],
      r[ 2 ] / r[ 0 ]
    };

    // Computing the first and the last plane
    int iMin = 0, iMax = 0;
    if( r[ 0 ] > 0.0 )
    {
      iMin = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alphaMin * r[ 0 ] - event1[ 0 ] ) / mp_sizeVox[ 0 ] );
      iMax = ::floor( m_toleranceX + ( event1[ 0 ] + alphaMax * r[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }
    else if( r[ 0 ] < 0.0 )
    {
      iMin = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alphaMax * r[ 0 ] - event1[ 0 ] ) / mp_sizeVox[ 0 ] );
      iMax = ::floor( m_toleranceX + ( event1[ 0 ] + alphaMin * r[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = -mp_halfFOV[ 0 ] + ( mp_sizeVox[ 0 ] * 0.5 ) - event1[ 0 ];

    //// tof : erf of previous voxel edge - Gaussian center
    //HPFLTNB previous_edge_erf =  erf( ( ( -mp_halfFOV[ 0 ] + iMin * mp_sizeVox[ 0 ] - event1[ 0 ] ) * factor_for_tof - lor_tof_center) / tof_sigma_sqrt2 );
    //HPFLTNB next_edge_erf = 0.;
    HPFLTNB tof_weight = 0.;
    
    // Loop over all the crossing planes
    for( int i = iMin; i < iMax; ++i )
    {
      // Computing position on crossed plane in term of grid spacing
      HPFLTNB const step = offset + i * mp_sizeVox[ 0 ];
      pos[ 1 ] = event1[ 1 ] + step * ri[ 1 ];
      pos[ 2 ] = event1[ 2 ] + step * ri[ 2 ];

      // Find the index in the image
      index[ 0 ] = ::floor( m_toleranceY + ( pos[ 1 ] - m_boundY ) / mp_sizeVox[ 1 ] );
      index[ 1 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - m_boundZ ) / mp_sizeVox[ 2 ] );

      // Bilinear interpolation component (using floor)
      wfl[ 0 ] = pos[ 1 ] - ( m_boundY + index[ 0 ] * mp_sizeVox[ 1 ] );
      wfl[ 1 ] = pos[ 2 ] - ( m_boundZ + index[ 1 ] * mp_sizeVox[ 2 ] );
      wfl[ 0 ] /= mp_sizeVox[ 1 ];
      wfl[ 1 ] /= mp_sizeVox[ 2 ];

      // Bilinear interpolation component (using ceil)
      wcl[ 0 ] = 1.0 - wfl[ 0 ];
      wcl[ 1 ] = 1.0 - wfl[ 1 ];

      // Final coefficients
      w[ 0 ] = wcl[ 0 ] * wcl[ 1 ];
      w[ 1 ] = wcl[ 0 ] * wfl[ 1 ];
      w[ 2 ] = wfl[ 0 ] * wfl[ 1 ];
      w[ 3 ] = wfl[ 0 ] * wcl[ 1 ];

      // Check Y bounds
      if( index[ 0 ] <= 0 ) limitY1 = 0;
      if( index[ 0 ] >= mp_nbVox[ 1 ] ) limitY2 = 0;

      // Check Z bounds
      if( index[ 1 ] <= 0  ) limitZ1 = 0;
      if( index[ 1 ] >= mp_nbVox[ 2 ] ) limitZ2 = 0;

      // Storing values and indices
      finalIdx = i + ( index[ 0 ] - 1 ) * mp_nbVox[ 0 ] + ( index[ 1 ] - 1 ) * m_nbVoxXY;

      tof_weight = 0.;
      if (m_TOFWeightingFcnPrecomputedFlag)
      {
        // fetch the value of the precomputed TOF weighting function (centered at the TOF measurement)
        // nearest to the projection of the voxel center onto the LOR
        INTNB temp = m_TOFWeightingFcnNbSamples/2 + (INTNB)round( ( step * factor_for_tof - lor_tof_center ) * m_TOFPrecomputedSamplingFactor );
        if (temp>=0 && temp<m_TOFWeightingFcnNbSamples) tof_weight = mp_TOFWeightingFcn[temp];
      }
      else
      {
        //// tof : erf of next voxel edge - Gaussian center
        //next_edge_erf = erf( ( ( step +  mp_sizeVox[ 0 ] * 0.5 ) * factor_for_tof - lor_tof_center ) / tof_sigma_sqrt2 );
        //// integration of the Gaussian is done on the LOR portion matching the whole current voxel along X axis
        //tof_weight = 0.5 * ::fabs(previous_edge_erf - next_edge_erf) ;
        //// keeping record of the previous edge, so as to save 1 erf computation
        //previous_edge_erf = next_edge_erf;

        if (m_TOFBinProperProcessingFlag)
        {
          // integration between the edges of the current quantization TOF bin (centered at the TOF measurement)
          // of the Gaussian centered at the projection of the voxel center onto the LOR
          // normalized by the size of the bin so that the integral of the final TOF weighting function remains 1
          HPFLTNB temp_erf = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center - m_TOFBinSizeInMm/2. - step * factor_for_tof));
          HPFLTNB prev_erf = erf(temp_erf/tof_sigma_sqrt2);
          temp_erf = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center + m_TOFBinSizeInMm/2. - step * factor_for_tof));
          HPFLTNB new_erf = erf(temp_erf/tof_sigma_sqrt2);
          tof_weight  = 0.5 * fabs(new_erf - prev_erf) / m_TOFBinSizeInMm;
        }
        else
        {
          // value of the normalized Gaussian (centered at the TOF measurement)
          // at the projection of the voxel center onto the LOR
          HPFLTNB temp = ( step * factor_for_tof - lor_tof_center ) / tof_sigma;
          tof_weight =  exp(- 0.5 * temp * temp ) * m_TOFGaussianNormCoef;
        }
      }

      // check if the main voxel is relevant
      if (!m_hasMask || (limitY1 && limitZ1 && index[ 0 ] <= mp_nbVox[ 1 ] && index[ 1 ] <= mp_nbVox[ 2 ] && mp_mask[finalIdx]))
      {
        ap_ProjectionLine->AddVoxel(a_direction, finalIdx * limitY1 * limitZ1, w[ 0 ] * weight * tof_weight * limitY1 * limitZ1);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + m_nbVoxXY ) * limitY1 * limitZ2, w[ 1 ] * weight * tof_weight * limitY1 * limitZ2);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + mp_nbVox[ 0 ] + m_nbVoxXY ) * limitY2 * limitZ2, w[ 2 ] * weight * tof_weight * limitY2 * limitZ2);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + mp_nbVox[ 0 ] ) * limitY2 * limitZ1, w[ 3 ] * weight * tof_weight * limitY2 * limitZ1);
      }

      limitY1 = 1.0; limitY2 = 1.0;
      limitZ1 = 1.0; limitZ2 = 1.0;
    }
  }
  else
  {
    // Computing the parametric values
    // For the X-axis
    // We introduce 1 voxel size around Y-axis
    if( r[ 0 ] != 0.0 )
    {
      HPFLTNB const alphaX_0 = ( tof_edge_low[ 0 ] - mp_sizeVox[ 0 ] * 0.5 - event1[ 0 ] )
        / r[ 0 ];
      HPFLTNB const alphaX_1 = ( tof_edge_high[ 0 ] + mp_sizeVox[ 0 ] * 0.5 - event1[ 0 ] )
        / r[ 0 ];
      alphaMin = std::max( alphaMin, std::min( alphaX_0, alphaX_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaX_0, alphaX_1 ) );
    }

    // For the Y-axis
    // Make a shift on first and last plane (like an offset)
    // We stay in the limit of the image space
    HPFLTNB const alphaY_0 = ( tof_edge_low[ 1 ] - event1[ 1 ] ) / r[ 1 ];
    HPFLTNB const alphaY_1 = ( tof_edge_high[ 1 ] - event1[ 1 ] ) / r[ 1 ];
    alphaMin = std::max( alphaMin, std::min( alphaY_0, alphaY_1 ) );
    alphaMax = std::min( alphaMax, std::max( alphaY_0, alphaY_1 ) );

    // For the Z-axis 
    // We introduce 1 voxel size around Z-axis
    // Make a shift on first and last plane (like an offset)
    if( r[ 2 ] != 0.0 )
    {
      HPFLTNB const alphaZ_0 = ( tof_edge_low[ 2 ] - mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      HPFLTNB const alphaZ_1 = ( tof_edge_high[ 2 ] + mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      alphaMin = std::max( alphaMin, std::min( alphaZ_0, alphaZ_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaZ_0, alphaZ_1 ) );
    }

    if( alphaMax <= alphaMin ) return 0;

    if( r[ 0 ] == 0.0 ) if( event1[ 0 ] < -mp_halfFOV[ 0 ] || event1[ 0 ] > mp_halfFOV[ 0 ] ) return 0;
    if( r[ 2 ] == 0.0 ) if( event1[ 2 ] < -mp_halfFOV[ 2 ] || event1[ 2 ] > mp_halfFOV[ 2 ] ) return 0;

    // Computing weight of normalization
    HPFLTNB const factor( ::fabs( r[ 1 ] ) / lor_length );
    HPFLTNB const factor_for_tof( lor_length / r[ 1 ] );
    HPFLTNB const weight( mp_sizeVox[ 1 ] / factor );

    // Computing the increment
    HPFLTNB const ri[ 3 ] = {
      r[ 0 ] / r[ 1 ],
      1.0, // r[ 1 ] / r[ 1 ]
      r[ 2 ] / r[ 1 ]
    };

    // Computing the first and the last plane
    int jMin = 0, jMax = 0;
    if( r[ 1 ] > 0.0 )
    {
      jMin = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alphaMin * r[ 1 ] - event1[ 1 ] ) / mp_sizeVox[ 1 ] );
      jMax = ::floor( m_toleranceY + ( event1[ 1 ] + alphaMax * r[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }
    else if( r[ 1 ] < 0.0 )
    {
      jMin = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alphaMax * r[ 1 ] - event1[ 1 ] ) / mp_sizeVox[ 1 ] );
      jMax = ::floor( m_toleranceY + ( event1[ 1 ] + alphaMin * r[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }

    //// tof : erf of previous voxel edge - Gaussian center
    //HPFLTNB previous_edge_erf =  erf( ( ( -mp_halfFOV[ 1 ] + jMin * mp_sizeVox[ 1 ] - event1[ 1 ]) * factor_for_tof - lor_tof_center) / tof_sigma_sqrt2 );
    //HPFLTNB next_edge_erf = 0.;
    HPFLTNB tof_weight = 0.;

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = -mp_halfFOV[ 1 ] + ( mp_sizeVox[ 1 ] * 0.5 ) - event1[ 1 ];

    // Loop over all the crossing planes
    for( int j = jMin; j < jMax; ++j )
    {
      // Computing position on crossed plane in term of grid spacing
      HPFLTNB const step = offset + j * mp_sizeVox[ 1 ];
      pos[ 0 ] = event1[ 0 ] + step * ri[ 0 ];
      pos[ 2 ] = event1[ 2 ] + step * ri[ 2 ];

      // Find the index in the image
      index[ 0 ] = ::floor( m_toleranceX + ( pos[ 0 ] - m_boundX ) / mp_sizeVox[ 0 ] );
      index[ 1 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - m_boundZ ) / mp_sizeVox[ 2 ] );

      // Bilinear interpolation component (using floor)
      wfl[ 0 ] = pos[ 0 ] - ( m_boundX + index[ 0 ] * mp_sizeVox[ 0 ] );
      wfl[ 1 ] = pos[ 2 ] - ( m_boundZ + index[ 1 ] * mp_sizeVox[ 2 ] );
      wfl[ 0 ] /= mp_sizeVox[ 0 ];
      wfl[ 1 ] /= mp_sizeVox[ 2 ];

      // Bilinear interpolation component (using ceil)
      wcl[ 0 ] = 1.0 - wfl[ 0 ];
      wcl[ 1 ] = 1.0 - wfl[ 1 ];

      // Final coefficients
      w[ 0 ] = wcl[ 0 ] * wcl[ 1 ];
      w[ 1 ] = wcl[ 0 ] * wfl[ 1 ];
      w[ 2 ] = wfl[ 0 ] * wfl[ 1 ];
      w[ 3 ] = wfl[ 0 ] * wcl[ 1 ];

      // Check X bounds
      if( index[ 0 ] <= 0 ) limitX1 = 0;
      if( index[ 0 ] >= mp_nbVox[ 0 ] ) limitX2 = 0;

      // Check Z bounds
      if( index[ 1 ] <= 0 ) limitZ1 = 0;
      if( index[ 1 ] >= mp_nbVox[ 2 ] ) limitZ2 = 0;

      // Storing values and indices
      finalIdx = ( index[ 0 ] - 1 ) + j * mp_nbVox[ 0 ] + ( index[ 1 ] - 1 ) * m_nbVoxXY;

      tof_weight = 0.;
      if (m_TOFWeightingFcnPrecomputedFlag)
      {
        // fetch the value of the precomputed TOF weighting function (centered at the TOF measurement)
        // nearest to the projection of the voxel center onto the LOR
        INTNB temp = m_TOFWeightingFcnNbSamples/2 + (INTNB)round( ( step * factor_for_tof - lor_tof_center) * m_TOFPrecomputedSamplingFactor );
        if (temp>=0 && temp<m_TOFWeightingFcnNbSamples) tof_weight = mp_TOFWeightingFcn[temp];
      }
      else
      {
        //// tof : erf of next voxel edge - Gaussian center
        //next_edge_erf = erf( ( ( step +  mp_sizeVox[ 1 ] * 0.5 ) * factor_for_tof - lor_tof_center ) / tof_sigma_sqrt2 );
        //// integration of the Gaussian is done on the LOR portion matching the whole current voxel along X axis
        //tof_weight = 0.5 * ::fabs(previous_edge_erf - next_edge_erf) ;
        //// keeping record of the previous edge, so as to save 1 erf computation
        //previous_edge_erf = next_edge_erf;

        if (m_TOFBinProperProcessingFlag)
        {
          // integration between the edges of the current quantization TOF bin (centered at the TOF measurement)
          // of the Gaussian centered at the projection of the voxel center onto the LOR
          // normalized by the size of the bin so that the integral of the final TOF weighting function remains 1
          HPFLTNB temp_erf = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center - m_TOFBinSizeInMm/2. - step * factor_for_tof));
          HPFLTNB prev_erf = erf(temp_erf/tof_sigma_sqrt2);
          temp_erf = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center + m_TOFBinSizeInMm/2. - step * factor_for_tof));
          HPFLTNB new_erf = erf(temp_erf/tof_sigma_sqrt2);
          tof_weight  = 0.5 * fabs(new_erf - prev_erf) / m_TOFBinSizeInMm;
        }
        else
        {
          // value of the normalized Gaussian (centered at the TOF measurement)
          // at the projection of the voxel center onto the LOR
          HPFLTNB temp = ( step * factor_for_tof - lor_tof_center ) / tof_sigma;
          tof_weight =  exp(- 0.5 * temp * temp ) * m_TOFGaussianNormCoef;
        }
      }

      // check if the main voxel is relevant
      if (!m_hasMask || (limitX1 && limitZ1 && index[ 0 ] <= mp_nbVox[ 0 ] &&  index[ 1 ] <= mp_nbVox[ 2 ] && mp_mask[finalIdx]))
      {
        ap_ProjectionLine->AddVoxel(a_direction, finalIdx * limitX1 * limitZ1, w[ 0 ] * weight * tof_weight * limitX1 * limitZ1);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + m_nbVoxXY ) * limitX1 * limitZ2, w[ 1 ] * weight * tof_weight * limitX1 * limitZ2);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + 1 + m_nbVoxXY ) * limitX2 * limitZ2, w[ 2 ] * weight * tof_weight * limitX2 * limitZ2);
        ap_ProjectionLine->AddVoxel(a_direction, ( finalIdx + 1 ) * limitX2 * limitZ1, w[ 3 ] * weight * tof_weight * limitX2 * limitZ1);
      }

      limitX1 = 1; limitX2 = 1;
      limitZ1 = 1; limitZ2 = 1;
    }
  }
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorJoseph::ProjectTOFHistogram(int a_direction, oProjectionLine* ap_ProjectionLine)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorJoseph::ProjectTOFHistogram() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorJoseph::Project with TOF bins -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // Get event positions
  FLTNB* event1 = ap_ProjectionLine->GetPosition1();
  FLTNB* event2 = ap_ProjectionLine->GetPosition2();

  // Distance between point event1 and event2
  HPFLTNB const r[ 3 ] = {
    event2[ 0 ] - event1[ 0 ],
    event2[ 1 ] - event1[ 1 ],
    event2[ 2 ] - event1[ 2 ]
  };

  // LOR length
  HPFLTNB lor_length = ap_ProjectionLine->GetLength();
  HPFLTNB lor_length_half = lor_length * 0.5;
  
  // Get TOF info
  INTNB tof_nb_bins = ap_ProjectionLine->GetNbTOFBins();
  INTNB tof_half_nb_bins = tof_nb_bins/2;

  // TOF Gaussian standard deviation and truncation
  HPFLTNB tof_sigma = m_TOFResolutionInMm / TWO_SQRT_TWO_LN_2;
  HPFLTNB tof_sigma_sqrt2 = sqrt(2.)*tof_sigma;
  HPFLTNB tof_half_span = 0.;
  if (m_TOFWeightingFcnPrecomputedFlag) tof_half_span = ((HPFLTNB)m_TOFWeightingFcnNbSamples)/(2.*m_TOFPrecomputedSamplingFactor);
  else tof_half_span = tof_sigma * m_TOFNbSigmas;

  // if integration of the Gaussian over the TOF bin, need for help variables to save calls to erf
  HPFLTNB prev_erf = 0., new_erf = 0.;

  // minimum and maximum TOF bins, the whole span
  INTNB tof_bin_last = tof_half_nb_bins;
  INTNB tof_bin_first = -tof_half_nb_bins;

  // distance between the first event1 and the center of a TOF bin along the LOR
  HPFLTNB lor_tof_center = 0.;
  // the sum of all TOF bin weights for a voxel
  //HPFLTNB tof_norm_coef = 0.;

  // the first and the last relevant TOF bins for a voxel
  INTNB tof_bin_first_for_voxel = 0, tof_bin_last_for_voxel = 0;

/*
  // Square of r
  HPFLTNB const r2[ 3 ] = {
    r[ 0 ] * r[ 0 ],
    r[ 1 ] * r[ 1 ],
    r[ 2 ] * r[ 2 ]
  };
*/

  // Find the first and last intersecting plane in X using the parametric
  // values alphaMin and alphaMax
  HPFLTNB alphaMin = 0.0, alphaMax = 1.0;

  // Variables for Joseph
  HPFLTNB pos[ 3 ] = { 0.0, 0.0, 0.0 }; // Position of point in image space
  HPFLTNB wfl[ 2 ] = { 0.0, 0.0 }; // Weight floor
  HPFLTNB wcl[ 2 ] = { 0.0, 0.0 }; // Weight ceil
  HPFLTNB w[ 4 ] = { 0.0, 0.0, 0.0, 0.0 }; // Interpolation weight
  int16_t index[ 2 ] = { 0, 0 }; // index in the image space
  int32_t finalIdx = 0; // final index
  int8_t limitX1 = 1; int8_t limitX2 = 1;
  int8_t limitY1 = 1; int8_t limitY2 = 1;
  int8_t limitZ1 = 1; int8_t limitZ2 = 1;

  // Condition on the largest component of r
  if( ::fabs( r[ 0 ] ) > ::fabs( r[ 1 ] ) )
  {
    // Computing the parametric values
    // For the X-axis
    // We stay in the limit of the image space
    HPFLTNB const alphaX_0 = ( -mp_halfFOV[ 0 ] - event1[ 0 ] ) / r[ 0 ];
    HPFLTNB const alphaX_1 = (  mp_halfFOV[ 0 ] - event1[ 0 ] ) / r[ 0 ];
    alphaMin = std::max( alphaMin, std::min( alphaX_0, alphaX_1 ) );
    alphaMax = std::min( alphaMax, std::max( alphaX_0, alphaX_1 ) );

    // For the Y-axis
    // We introduce 1 voxel size around Y-axis
    // Make a shift on first and last plane (like an offset)
    if( r[ 1 ] != 0.0 )
    {
      HPFLTNB const alphaY_0 = ( -mp_halfFOV[ 1 ] - mp_sizeVox[ 1 ] * 0.5 - event1[ 1 ] )
        / r[ 1 ];
      HPFLTNB const alphaY_1 = ( mp_halfFOV[ 1 ] + mp_sizeVox[ 1 ] * 0.5 - event1[ 1 ] )
        / r[ 1 ];
      alphaMin = std::max( alphaMin, std::min( alphaY_0, alphaY_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaY_0, alphaY_1 ) );
    }

    // For the Z-axis 
    // We introduce 1 voxel size around Z-axis
    // Make a shift on first and last plane (like an offset)
    if( r[ 2 ] != 0.0 )
    {
      HPFLTNB const alphaZ_0 = ( -mp_halfFOV[ 2 ] - mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      HPFLTNB const alphaZ_1 = ( mp_halfFOV[ 2 ] + mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      alphaMin = std::max( alphaMin, std::min( alphaZ_0, alphaZ_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaZ_0, alphaZ_1 ) );
    }

    if( alphaMax <= alphaMin ) return 0;

    if( r[ 1 ] == 0.0 ) if( event1[ 1 ] < -mp_halfFOV[ 1 ] || event1[ 1 ] > mp_halfFOV[ 1 ] ) return 0;
    if( r[ 2 ] == 0.0 ) if( event1[ 2 ] < -mp_halfFOV[ 2 ] || event1[ 2 ] > mp_halfFOV[ 2 ] ) return 0;

    // temporary storage for TOF bin weights
    // allocation after potential returns
    HPFLTNB* tof_weights_temp = new HPFLTNB[tof_nb_bins];
    for (INTNB tof_bin=0; tof_bin<tof_nb_bins; tof_bin++) tof_weights_temp[tof_bin] = 0.;

    // Computing weight of normalization
    HPFLTNB const factor( ::fabs( r[ 0 ] ) / lor_length );
    HPFLTNB const factor_for_tof( lor_length / r[ 0 ] );
    HPFLTNB const weight( mp_sizeVox[ 0 ] / factor );

    // Computing the increment
    HPFLTNB const ri[ 3 ] = {
      1.0, // r[ 0 ] / r[ 0 ]
      r[ 1 ] / r[ 0 ],
      r[ 2 ] / r[ 0 ]
    };

    // Computing the first and the last plane
    int iMin = 0, iMax = 0;
    if( r[ 0 ] > 0.0 )
    {
      iMin = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alphaMin * r[ 0 ] - event1[ 0 ] ) / mp_sizeVox[ 0 ] );
      iMax = ::floor( m_toleranceX + ( event1[ 0 ] + alphaMax * r[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }
    else if( r[ 0 ] < 0.0 )
    {
      iMin = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alphaMax * r[ 0 ] - event1[ 0 ] ) / mp_sizeVox[ 0 ] );
      iMax = ::floor( m_toleranceX + ( event1[ 0 ] + alphaMin * r[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = -mp_halfFOV[ 0 ] + ( mp_sizeVox[ 0 ] * 0.5 ) - event1[ 0 ];

    // Loop over all the crossing planes
    for( int i = iMin; i < iMax; ++i )
    {
      // Computing position on crossed plane in term of grid spacing
      HPFLTNB const step = offset + i * mp_sizeVox[ 0 ];
      pos[ 1 ] = event1[ 1 ] + step * ri[ 1 ];
      pos[ 2 ] = event1[ 2 ] + step * ri[ 2 ];

      // Find the index in the image
      index[ 0 ] = ::floor( m_toleranceY + ( pos[ 1 ] - m_boundY ) / mp_sizeVox[ 1 ] );
      index[ 1 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - m_boundZ ) / mp_sizeVox[ 2 ] );

      finalIdx = i + ( index[ 0 ] - 1 ) * mp_nbVox[ 0 ] + ( index[ 1 ] - 1 ) * m_nbVoxXY;

      // check if the main voxel is relevant
      if (m_hasMask && index[ 0 ] > 0 && index[ 1 ] > 0 && index[ 0 ] <= mp_nbVox[ 1 ] &&  index[ 1 ] <= mp_nbVox[ 2 ] && !mp_mask[finalIdx]) continue;

      // Compute the first and the last relevant TOF bin for this voxel
      if (!m_TOFWeightingFcnPrecomputedFlag && m_TOFBinProperProcessingFlag)
      {
        // taking the TOF bins reached by the truncated Gaussian centered at the voxel center projection
        tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(floor( (step * factor_for_tof - tof_half_span - lor_length_half - m_TOFBinSizeInMm/2.) / m_TOFBinSizeInMm )));
        tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(ceil( (step * factor_for_tof + tof_half_span - lor_length_half + m_TOFBinSizeInMm/2.) / m_TOFBinSizeInMm )));
      }
      else
      {
        // taking the TOF bins whose TOF weighting function, centered at the bin center, reaches the voxel center projection
        tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(ceil( (step * factor_for_tof - tof_half_span - lor_length_half ) / m_TOFBinSizeInMm )));
        tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(floor( (step * factor_for_tof + tof_half_span - lor_length_half) / m_TOFBinSizeInMm )));
      }

      lor_tof_center = lor_length_half + tof_bin_first_for_voxel * m_TOFBinSizeInMm;

      // shift tof bin indices from -/+ to 0:nbBins range
      tof_bin_first_for_voxel += tof_half_nb_bins;
      tof_bin_last_for_voxel += tof_half_nb_bins;

      // initialization of help variables for reducing calls to erf
      if (!m_TOFWeightingFcnPrecomputedFlag && m_TOFBinProperProcessingFlag)
      {
        // bound the integration to the Gaussian truncation
        HPFLTNB temp = std::min(tof_half_span,std::max(-tof_half_span, lor_tof_center - m_TOFBinSizeInMm/2. - step * factor_for_tof));
        prev_erf = erf(temp/tof_sigma_sqrt2);
      }

      // compute the TOF bin weights for the current voxel for all relevant TOF bins
      //tof_norm_coef = 0.;
      for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++)
      {
        if (m_TOFWeightingFcnPrecomputedFlag)
        {
          // fetch the value of the precomputed TOF weighting function (centered at the TOF bin center)
          // nearest to the projection of the voxel center onto the LOR
          INTNB temp = m_TOFWeightingFcnNbSamples/2 + (INTNB)round( (step * factor_for_tof - lor_tof_center) * m_TOFPrecomputedSamplingFactor);
          if (temp>=0 && temp<m_TOFWeightingFcnNbSamples) tof_weights_temp[tof_bin] = mp_TOFWeightingFcn[temp];
          // add the weight to the sum
          //tof_norm_coef += tof_weights_temp[tof_bin];
          // update TOF bin center along the LOR for the next TOF bin
          lor_tof_center += m_TOFBinSizeInMm;
        }
        else
        {
          if (m_TOFBinProperProcessingFlag)
          {
            // integration between the edges of the current TOF bin of the Gaussian centered at the projection of the voxel center onto the LOR
            // reuse of integration from the previous bin to save calls to erf
            // bound the integration to the Gaussian truncation
            HPFLTNB temp = std::min(tof_half_span,std::max(-tof_half_span, lor_tof_center + m_TOFBinSizeInMm/2. - step * factor_for_tof));
            new_erf = erf(temp/tof_sigma_sqrt2);
            tof_weights_temp[tof_bin] = 0.5 * fabs(new_erf - prev_erf);
            // add the weight to the sum
            //tof_norm_coef += tof_weights_temp[tof_bin];
            prev_erf = new_erf;
            // update TOF bin center along the LOR for the next TOF bin
            lor_tof_center += m_TOFBinSizeInMm;
          }
          else
          {
            // TOF weight = TOF bin size * value of the normalized Gaussian (centered at the TOF bin center) at the projection of the voxel center onto the LOR
            HPFLTNB temp = ( step * factor_for_tof - lor_tof_center) / tof_sigma;
            // save the weight temporarily
            tof_weights_temp[tof_bin] = exp(- 0.5 * temp * temp ) * m_TOFGaussianNormCoef * m_TOFBinSizeInMm;
            // add the weight to the sum
            //tof_norm_coef += tof_weights_temp[tof_bin];
            // update TOF bin center along the LOR for the next TOF bin
            lor_tof_center += m_TOFBinSizeInMm;
          }
        }
      }

      // Bilinear interpolation component (using floor)
      wfl[ 0 ] = pos[ 1 ] - ( m_boundY + index[ 0 ] * mp_sizeVox[ 1 ] );
      wfl[ 1 ] = pos[ 2 ] - ( m_boundZ + index[ 1 ] * mp_sizeVox[ 2 ] );
      wfl[ 0 ] /= mp_sizeVox[ 1 ];
      wfl[ 1 ] /= mp_sizeVox[ 2 ];

      // Bilinear interpolation component (using ceil)
      wcl[ 0 ] = 1.0 - wfl[ 0 ];
      wcl[ 1 ] = 1.0 - wfl[ 1 ];

      // Final coefficients
      w[ 0 ] = wcl[ 0 ] * wcl[ 1 ];
      w[ 1 ] = wcl[ 0 ] * wfl[ 1 ];
      w[ 2 ] = wfl[ 0 ] * wfl[ 1 ];
      w[ 3 ] = wfl[ 0 ] * wcl[ 1 ];

      // Check Y bounds
      if( index[ 0 ] <= 0 ) limitY1 = 0;
      if( index[ 0 ] >= mp_nbVox[ 1 ] ) limitY2 = 0;

      // Check Z bounds
      if( index[ 1 ] <= 0  ) limitZ1 = 0;
      if( index[ 1 ] >= mp_nbVox[ 2 ] ) limitZ2 = 0;
/*
      if (tof_norm_coef>0.)
      {
        // first normalize TOF bin weights so that they sum to 1
        for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++) tof_weights_temp[tof_bin] /= tof_norm_coef;
      }
*/
      // compute and write the final TOF bin projection coefficients for current voxels
      ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, finalIdx * limitY1 * limitZ1, w[ 0 ] * weight * limitY1 * limitZ1, tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
      ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, ( finalIdx + m_nbVoxXY ) * limitY1 * limitZ2, w[ 1 ] * weight * limitY1 * limitZ2 , tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
      ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, ( finalIdx + mp_nbVox[ 0 ] + m_nbVoxXY ) * limitY2 * limitZ2, w[ 2 ] * weight * limitY2 * limitZ2 , tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
      ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, ( finalIdx + mp_nbVox[ 0 ] ) * limitY2 * limitZ1, w[ 3 ] * weight * limitY2 * limitZ1 , tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
      

      limitY1 = 1.0; limitY2 = 1.0;
      limitZ1 = 1.0; limitZ2 = 1.0;
    }
    delete[] tof_weights_temp;
  }
  else
  {
    // Computing the parametric values
    // For the X-axis
    // We introduce 1 voxel size around Y-axis
    if( r[ 0 ] != 0.0 )
    {
      HPFLTNB const alphaX_0 = ( -mp_halfFOV[ 0 ] - mp_sizeVox[ 0 ] * 0.5 - event1[ 0 ] )
        / r[ 0 ];
      HPFLTNB const alphaX_1 = ( mp_halfFOV[ 0 ] + mp_sizeVox[ 0 ] * 0.5 - event1[ 0 ] )
        / r[ 0 ];
      alphaMin = std::max( alphaMin, std::min( alphaX_0, alphaX_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaX_0, alphaX_1 ) );
    }

    // For the Y-axis
    // Make a shift on first and last plane (like an offset)
    // We stay in the limit of the image space
    HPFLTNB const alphaY_0 = ( -mp_halfFOV[ 1 ] - event1[ 1 ] ) / r[ 1 ];
    HPFLTNB const alphaY_1 = ( mp_halfFOV[ 1 ] - event1[ 1 ] ) / r[ 1 ];
    alphaMin = std::max( alphaMin, std::min( alphaY_0, alphaY_1 ) );
    alphaMax = std::min( alphaMax, std::max( alphaY_0, alphaY_1 ) );

    // For the Z-axis 
    // We introduce 1 voxel size around Z-axis
    // Make a shift on first and last plane (like an offset)
    if( r[ 2 ] != 0.0 )
    {
      HPFLTNB const alphaZ_0 = ( -mp_halfFOV[ 2 ] - mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      HPFLTNB const alphaZ_1 = ( mp_halfFOV[ 2 ] + mp_sizeVox[ 2 ] * 0.5 - event1[ 2 ] )
        / r[ 2 ];
      alphaMin = std::max( alphaMin, std::min( alphaZ_0, alphaZ_1 ) );
      alphaMax = std::min( alphaMax, std::max( alphaZ_0, alphaZ_1 ) );
    }

    if( alphaMax <= alphaMin ) return 0;

    if( r[ 0 ] == 0.0 ) if( event1[ 0 ] < -mp_halfFOV[ 0 ] || event1[ 0 ] > mp_halfFOV[ 0 ] ) return 0;
    if( r[ 2 ] == 0.0 ) if( event1[ 2 ] < -mp_halfFOV[ 2 ] || event1[ 2 ] > mp_halfFOV[ 2 ] ) return 0;

    // temporary storage for TOF bins Gaussian integrals over the current voxel, allocation after potential returns
    HPFLTNB* tof_weights_temp = new HPFLTNB[tof_nb_bins];
    for (INTNB tof_bin=0; tof_bin<tof_nb_bins; tof_bin++) tof_weights_temp[tof_bin] = 0.;
    
    // Computing weight of normalization
    HPFLTNB const factor( ::fabs( r[ 1 ] ) / lor_length );
    HPFLTNB const factor_for_tof( lor_length / r[ 1 ] );
    HPFLTNB const weight( mp_sizeVox[ 1 ] / factor );

    // Computing the increment
    HPFLTNB const ri[ 3 ] = {
      r[ 0 ] / r[ 1 ],
      1.0, // r[ 1 ] / r[ 1 ]
      r[ 2 ] / r[ 1 ]
    };

    // Computing the first and the last plane
    int jMin = 0, jMax = 0;
    if( r[ 1 ] > 0.0 )
    {
      jMin = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alphaMin * r[ 1 ] - event1[ 1 ] ) / mp_sizeVox[ 1 ] );
      jMax = ::floor( m_toleranceY + ( event1[ 1 ] + alphaMax * r[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }
    else if( r[ 1 ] < 0.0 )
    {
      jMin = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alphaMax * r[ 1 ] - event1[ 1 ] ) / mp_sizeVox[ 1 ] );
      jMax = ::floor( m_toleranceY + ( event1[ 1 ] + alphaMin * r[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = -mp_halfFOV[ 1 ] + ( mp_sizeVox[ 1 ] * 0.5 ) - event1[ 1 ];

    // Loop over all the crossing planes
    for( int j = jMin; j < jMax; ++j )
    {
      // Computing position on crossed plane in term of grid spacing
      HPFLTNB const step = offset + j * mp_sizeVox[ 1 ];
      pos[ 0 ] = event1[ 0 ] + step * ri[ 0 ];
      pos[ 2 ] = event1[ 2 ] + step * ri[ 2 ];

      // Find the index in the image
      index[ 0 ] = ::floor( m_toleranceX + ( pos[ 0 ] - m_boundX ) / mp_sizeVox[ 0 ] );
      index[ 1 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - m_boundZ ) / mp_sizeVox[ 2 ] );
  
      finalIdx = ( index[ 0 ] - 1 ) + j * mp_nbVox[ 0 ] + ( index[ 1 ] - 1 ) * m_nbVoxXY;

      // check if the main voxel is relevant
      if (m_hasMask && (index[ 0 ] > 0 && index[ 1 ] > 0 && index[ 0 ] <= mp_nbVox[ 0 ] &&  index[ 1 ] <= mp_nbVox[ 2 ]) && !mp_mask[finalIdx]) continue;

      // Compute the first and the last relevant TOF bin for this voxel
      if (!m_TOFWeightingFcnPrecomputedFlag && m_TOFBinProperProcessingFlag)
      {
        // taking the TOF bins reached by the truncated Gaussian centered at the voxel center projection
        tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(floor( (step * factor_for_tof - tof_half_span - lor_length_half - m_TOFBinSizeInMm/2.) / m_TOFBinSizeInMm )));
        tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(ceil( (step * factor_for_tof + tof_half_span - lor_length_half + m_TOFBinSizeInMm/2.) / m_TOFBinSizeInMm )));
      }
      else
      {
        // taking the TOF bins whose TOF uncertainty function (centered at the bin center) reaches the voxel center projection 
        tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(ceil( (step * factor_for_tof - tof_half_span - lor_length_half ) / m_TOFBinSizeInMm )));
        tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(floor( (step * factor_for_tof + tof_half_span - lor_length_half) / m_TOFBinSizeInMm )));
      }

      lor_tof_center = lor_length_half + tof_bin_first_for_voxel * m_TOFBinSizeInMm;

      // if min/max TOF bins for this voxel do not fall into the total available TOF bins range, skip
      if(tof_bin_first_for_voxel>tof_bin_last || tof_bin_last_for_voxel<tof_bin_first) continue;
      
      // shift tof bin indices from -/+ to 0:nbBins range
      tof_bin_first_for_voxel += tof_half_nb_bins;
      tof_bin_last_for_voxel += tof_half_nb_bins;

      // initialization of help variables for reducing calls to erf
      if (!m_TOFWeightingFcnPrecomputedFlag && m_TOFBinProperProcessingFlag)
      {
        // bound the integration to the Gaussian truncation
        HPFLTNB temp = std::min(tof_half_span,std::max(-tof_half_span, lor_tof_center - m_TOFBinSizeInMm/2. - step * factor_for_tof));
        prev_erf = erf(temp/tof_sigma_sqrt2);
      }

      // compute TOF bin weights for the current voxel for all relevant TOF bins
      //tof_norm_coef = 0.;
      for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++)
      {
        if (m_TOFWeightingFcnPrecomputedFlag)
        {
          // fetch the value of the precomputed TOF weighting function (centered at the TOF bin center)
          // nearest to the projection of the voxel center onto the LOR
          INTNB temp = m_TOFWeightingFcnNbSamples/2 + (INTNB)round( (step * factor_for_tof - lor_tof_center) * m_TOFPrecomputedSamplingFactor );
          if (temp>=0 && temp<m_TOFWeightingFcnNbSamples) tof_weights_temp[tof_bin] = mp_TOFWeightingFcn[temp];
          // add the weight to the sum for normalization
          //tof_norm_coef += tof_weights_temp[tof_bin];
          // update TOF center along the LOR for the next TOF bin
          lor_tof_center += m_TOFBinSizeInMm;
        }
        else
        {
          if (m_TOFBinProperProcessingFlag)
          {
            // integration between the edges of the current TOF bin of the Gaussian centered at the projection of the voxel center onto the LOR
            // reuse of integration from the previous bin to save calls to erf
            // bound the integration to the Gaussian truncation
            HPFLTNB temp = std::min(tof_half_span,std::max(-tof_half_span, lor_tof_center + m_TOFBinSizeInMm/2. - step * factor_for_tof));
            new_erf = erf(temp/tof_sigma_sqrt2);
            tof_weights_temp[tof_bin] = 0.5 * fabs(new_erf - prev_erf);
            // add the weight to the sum for normalization
            //tof_norm_coef += tof_weights_temp[tof_bin];
            prev_erf = new_erf;
            // update TOF center along the LOR for the next TOF bin
            lor_tof_center += m_TOFBinSizeInMm;
          }
          else
          {
            // TOF weight = TOF bin size * value of the normalized Gaussian (centered at the TOF bin center) at the projection of the voxel center onto the LOR
            HPFLTNB temp = ( step * factor_for_tof - lor_tof_center) / tof_sigma;
            // save the weight temporarily
            tof_weights_temp[tof_bin] = exp(- 0.5 * temp * temp ) * m_TOFGaussianNormCoef * m_TOFBinSizeInMm;
            // add the weight to the sum
            //tof_norm_coef += tof_weights_temp[tof_bin];
            // update TOF bin center along the LOR for the next TOF bin
            lor_tof_center += m_TOFBinSizeInMm;
          }
        }
      }

      // Bilinear interpolation component (using floor)
      wfl[ 0 ] = pos[ 0 ] - ( m_boundX + index[ 0 ] * mp_sizeVox[ 0 ] );
      wfl[ 1 ] = pos[ 2 ] - ( m_boundZ + index[ 1 ] * mp_sizeVox[ 2 ] );
      wfl[ 0 ] /= mp_sizeVox[ 0 ];
      wfl[ 1 ] /= mp_sizeVox[ 2 ];

      // Bilinear interpolation component (using ceil)
      wcl[ 0 ] = 1.0 - wfl[ 0 ];
      wcl[ 1 ] = 1.0 - wfl[ 1 ];

      // Final coefficients
      w[ 0 ] = wcl[ 0 ] * wcl[ 1 ];
      w[ 1 ] = wcl[ 0 ] * wfl[ 1 ];
      w[ 2 ] = wfl[ 0 ] * wfl[ 1 ];
      w[ 3 ] = wfl[ 0 ] * wcl[ 1 ];

      // Check X bounds
      if( index[ 0 ] <= 0 ) limitX1 = 0;
      if( index[ 0 ] >= mp_nbVox[ 0 ] ) limitX2 = 0;

      // Check Z bounds
      if( index[ 1 ] <= 0 ) limitZ1 = 0;
      if( index[ 1 ] >= mp_nbVox[ 2 ] ) limitZ2 = 0;

      // Storing values and indices
/*
      if (tof_norm_coef>0.)
      {
        // first normalize TOF bin weights so that they sum to 1
        for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++) tof_weights_temp[tof_bin] /= tof_norm_coef;
      }
*/
      // compute and write the final TOF bin projection coefficients for current voxels
      ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, finalIdx * limitX1 * limitZ1, w[ 0 ] * weight * limitX1 * limitZ1, tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
      ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, ( finalIdx + m_nbVoxXY ) * limitX1 * limitZ2, w[ 1 ] * weight * limitX1 * limitZ2, tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
      ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, ( finalIdx + 1 + m_nbVoxXY ) * limitX2 * limitZ2, w[ 2 ] * weight * limitX2 * limitZ2, tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
      ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, ( finalIdx + 1 ) * limitX2 * limitZ1, w[ 3 ] * weight * limitX2 * limitZ1, tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);

      limitX1 = 1; limitX2 = 1;
      limitZ1 = 1; limitZ2 = 1;
    }
    delete[] tof_weights_temp;
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
