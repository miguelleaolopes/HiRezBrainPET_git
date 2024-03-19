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
  \ingroup datafile

  \brief Implementation of class iEventHistoCT
*/

#include "iEventHistoCT.hh"
#include "vDataFile.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventHistoCT::iEventHistoCT() : iEventCT()
{
  m_dataType = TYPE_CT;
  m_dataSpec = SPEC_TRANSMISSION;
  m_dataMode = MODE_HISTOGRAM;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventHistoCT::~iEventHistoCT() {}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventHistoCT::Describe()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("iEventHistoCT::Describe() -> Display contents" << endl);
  Cout("Time: " << m_timeInMs << " ms" << endl);
  Cout("Number of lines: " << m_nbLines << endl);
  for (uint16_t l=0; l<m_nbLines; l++) Cout("  --> ID1: " << mp_ID1[l] << " | ID2: " << mp_ID2[l] << endl);
  Cout("Bin value: " << m_eventValue << endl);
  Cout("Scatter rate: " << m_eventScatRate << endl);
  Cout("Blank value: " << m_eventBlankValue << endl);
  Cout(flush);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
