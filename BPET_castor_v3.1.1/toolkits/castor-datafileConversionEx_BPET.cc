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
  \ingroup utils
  \brief This program is an example code providing guidance regarding how to
         generate a PET CASToR datafile from any system datafile. \n
         It does not perform any conversion by itself and must be adjusted
         to the conversion of any system dataset. It uses a GATE model
         of a General Electric PET system as example. \n
         It implements the required calls to CASToR DataFile and Event
         objects and functions in order to write a dataset in the CASToR
         format, with or without optionnal correction factors (scatter,
         random, normalization coefficients, etc).
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "iDataFilePET.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"
#include <iomanip> //std::setprecision

#ifdef CASTOR_ROOT
#ifdef _WIN32
#include "Windows4Root.h"
#endif
#include "TROOT.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#endif

/*!
  \fn      ShowHelp()
  \brief   Display main command line options for castor-datafileConversionEx
*/
void ShowHelp(int a_returnCode)
{
  cout << endl;
  cout << "NOTE : This program is an example code providing guidance regarding how to generate a PET CASToR datafile from any system datafile " << endl;
  cout << "       It does NOT perform any conversion by itself and must be adjusted to the conversion of any system dataset !" << endl;
  cout << endl;
  cout << endl;
  cout << "Usage:  castor-datafileConversionEx -ih input_datafile : (if datafile is histogram mode)" << endl;
  cout << "                                    -il input_datafile : (if datafile is list-mode)" << endl;
  cout << "                                    -o  output_datafile" << endl;
  cout << "                                    -s  scanner_alias " << endl;
  cout << "                                   [-nc  norm_factors_file" << endl;
  cout << "                                   [-sc  scat_factors_file" << endl;
  cout << "                                   [-rc  rdm_factors_file" << endl;
  cout << "                                   [-ac  atn_factors_file" << endl;
  cout << "                                   [-cf  calibration factor" << endl;
  cout << "                                   [-ist isotope_alias" << endl;
  // cout << "                                [-poi  flag for POI]" << endl;
  // cout << "                                [-tof  flag for TOF]" << endl;
  cout << "                                   [-h    help" << endl;
  cout << "                                   [-vb   verbosity]" << endl;
  cout << endl;
  cout << endl;
  cout << "[Main settings]:" << endl;
  cout << "  -ih path_to_histo_datafile   : provide an input histogram datafile to convert" << endl;
  cout << "  -il path_to_listm_datafile   : provide an input list-mode datafile to convert" << endl;
  cout << "  -o  path_to_cstr_file.cdh    : provide the path to the output file will be created inside this folder (no default)" << endl;
  cout << "  -s  scanner_alias            : provide the name of the scanner used for to acquire the original data" << endl;
  cout << "                               : Must correspond to a .geom or .hscan file in the config/scanner repository." << endl;
  cout << endl;
  cout << "[Optional settings]:" << endl;
  cout << "  -nc norm_factors_file        : provide a file containing normalization correction factors" << endl;
  cout << "  -sc scat_factors_file        : provide a file containing scatter correction factors" << endl;
  cout << "  -rc rdm_factors_file         : provide a file containing random correction factors" << endl;
  cout << "  -ac atn_factors_file         : provide a file containing attenuation correction factors" << endl;
  cout << "  -cf calibration factor       : provide a calibration factor" << endl;
  cout << "  -ist isotope_alias           : provide alias of the isotope used in the input datafile" << endl;
  cout << "                                 Supported PET isotopes and their parameters are listed in config/misc/isotopes_pet" << endl;
  cout << "                                 Other isotopes can be added in the same file" << endl;
  cout << endl;
  cout << "[Miscellaneous settings]:" << endl;
  cout << "  -h                           : Print out this help page" << endl;
  cout << "  -vb                          : give the verbosity level, from 0 (no verbose) to above 5 (at the event level) (default: 1)." << endl;
  cout << "  -conf                        : Give the path to the CASToR config directory (default: located through the CASTOR_CONFIG environment variable)." << endl;
  cout << endl;
#ifdef CASTOR_VERSION
  cout << "  This program is part of the CASToR release version " << CASTOR_VERSION << "." << endl;
  cout << endl;
#endif
  Exit(a_returnCode);
}

/*
  Main program
*/
int main(int argc, char **argv)
{
// ============================================================================================================
// MPI stuff (we make all instances but the first one returning 0 directly)
// ============================================================================================================
#ifdef CASTOR_MPI
  int mpi_rank = 0;
  int mpi_size = 1;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_rank != 0)
    return 0;
#endif

  // No argument, then show help
  if (argc == 1)
    ShowHelp(EXIT_SUCCESS);

  // ============================================================================================================
  // Parameterized variables with their default values
  // ============================================================================================================

  // String containing the path to the input data filename provided by the user
  string path_to_input_file = "";

  // DataFile mode (histogram or list (default))
  int data_mode = MODE_LIST;

  // Path to user provided data
  string path_to_out_file = "";
  string path_to_data_file = "";
  string path_to_header_file = "";
  string path_to_norm_file = "";
  string path_to_scat_file = "";
  string path_to_rdm_file = "";
  string path_to_acf_file = "";
  FLTNB calib_factor = 1.;

  string scanner_alias = "";
  string istp_alias = "";

  // Verbosity
  int vb = 1;

  // Input datafile is histogram
  bool is_histogram_flag = 0;

  // Path to config directory
  string path_to_config_dir = "";

  // ============================================================================================================
  // Read command-line parameters
  // ============================================================================================================

  for (int i = 1; i < argc; i++)
  {
    string option = (string)argv[i];

    // Show help
    if (option == "-h")
      ShowHelp(EXIT_SUCCESS);
    
    // Provide an histogram datafile to convert
    if (option == "-ih")
    {
      if (!path_to_input_file.empty())
      {
        Cerr("***** castor-datafileConversionEx :: the following file names have already been provided (either -ih or -il option must be used): " << endl);
        Cout(path_to_input_file << endl);
        Exit(EXIT_FAILURE);
      }
      else
      {
        if (argv[i + 1] == NULL)
        {
          Cerr("***** castor-datafileConversionEx :: argument missing for option: " << option << endl);
          Exit(EXIT_FAILURE);
        }
        else
          path_to_input_file = (string)argv[i + 1];

        data_mode = MODE_HISTOGRAM;
        i++;
      }
    }

    // Provide a list-mode datafile to convert
    if (option == "-il")
    {
      if (!path_to_input_file.empty())
      {
        Cerr("***** castor-datafileConversionEx :: the following file names have already been provided (either -ih or -il option must be used): " << endl);
        Cout(path_to_input_file << endl);
        Exit(EXIT_FAILURE);
      }
      else
      {
        if (argv[i + 1] == NULL)
        {
          Cerr("***** castor-datafileConversionEx :: argument missing for option: " << option << endl);
          Exit(EXIT_FAILURE);
        }
        else
          path_to_input_file = (string)argv[i + 1];

        data_mode = MODE_LIST;
        i++;
      }
    }

    // Output CASToR datafile
    else if (option == "-o")
    {
      if (argv[i + 1] == NULL)
      {
        Cerr("***** castor-datafileConversionEx :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        path_to_out_file = (string)argv[i + 1];
      i++;
    }

    // Scanner alias
    else if (option == "-s")
    {
      if (argv[i + 1] == NULL)
      {
        Cerr("***** castor-datafileConversionEx :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        scanner_alias = (string)argv[i + 1];
      i++;
    }

    // Provide a norm correction file
    else if (option == "-nc")
    {
      if (argv[i + 1] == NULL)
      {
        Cerr("***** castor-datafileConversionEx :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        path_to_norm_file = (string)argv[i + 1];

      i++;
    }

    // Provide a scatter correction file
    else if (option == "-sc")
    {
      if (argv[i + 1] == NULL)
      {
        Cerr("***** castor-datafileConversionEx :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        path_to_scat_file = (string)argv[i + 1];

      i++;
    }

    // Provide a rdm correction file
    else if (option == "-rc")
    {
      if (argv[i + 1] == NULL)
      {
        Cerr("***** castor-datafileConversionEx :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        path_to_rdm_file = (string)argv[i + 1];

      i++;
    }

    // Provide an acf correction file
    else if (option == "-ac")
    {
      if (argv[i + 1] == NULL)
      {
        Cerr("***** castor-datafileConversionEx :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        path_to_acf_file = (string)argv[i + 1];

      i++;
    }

    // Provide a calibration factor
    else if (option == "-cf")
    {
      if (argv[i + 1] == NULL)
      {
        Cerr("***** castor-datafileConversionEx :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
      {
        string val_str = (string)argv[i + 1];
        if (ConvertFromString(val_str, &calib_factor))
        {
          Cerr("***** castor-datafileConversionEx :: Exception when trying to read'" << val_str << " for option: " << option << endl);
          Exit(EXIT_FAILURE);
        }
      }
      i++;
    }

    // Provide an isotope alias
    else if (option == "-ist")
    {
      if (argv[i + 1] == NULL)
      {
        Cerr("***** castor-datafileConversionEx :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        istp_alias = (string)argv[i + 1];

      i++;
    }

    // Verbosity level
    else if (option == "-vb")
    {
      if (i >= argc - 1)
      {
        Cerr("***** castor-datafileConversionEx :: Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      vb = atoi(argv[i + 1]);
      i++;
    }
    
    // Path to config directory
    else if (option=="-conf")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_config_dir = (string)argv[i+1];
      i++;
    }

    else
    {
      Cerr("***** castor-datafileConversionEx :: Unknown option '" << option << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // ============================================================================================================
  // Basic initialization checks (minimal initializations mandatory for the next steps)
  // ============================================================================================================

  // data files
  if (path_to_input_file.empty())
  {
    Cerr("***** castor-datafileConversionEx :: Please provide at least one data filename (-ih for histogram, -il for list-mode)" << endl);
    cout << "  -input filename   : give an input datafile (no default)." << endl;
    Exit(EXIT_FAILURE);
  }

  // output directory
  if (path_to_out_file.empty())
  {
    Cerr("***** castor-datafileConversionEx :: Please provide the output file name (-o)" << endl);
    Exit(EXIT_FAILURE);
  }

  // scanner
  if (scanner_alias.empty())
  {
    Cerr("***** castor-datafileConversionEx :: Please provide a scanner alias (-s) :" << endl);
    Cout("                                     It must correspond to a .geom or .hscan (user-made LUT) file in the config/scanner directory." << endl);
    ShowHelp(0);
    Exit(EXIT_FAILURE);
  }


  // ============================================================================================================
  // SOutputManager object initialisation:
  // (Used for log file)
  // ============================================================================================================

  sOutputManager *p_outputManager = sOutputManager::GetInstance();
  // Set verbose level
  p_outputManager->SetVerbose(vb);
  // Set MPI rank
  p_outputManager->SetMPIRank(0);

  // Set path to the config directory
  if (p_outputManager->CheckConfigDir(path_to_config_dir))
  {
    Cerr("***** castor-datafileConversionEx :: A problem occurred while checking for the config directory path !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize output directory and base name
  if (p_outputManager->InitOutputDirectory(path_to_out_file, ""))
  {
    Cerr("***** castor-datafileConversionEx :: A problem occurred while initializing output directory !" << endl);

    Exit(EXIT_FAILURE);
  }
  // Log command line
  if (p_outputManager->LogCommandLine(argc, argv))
  {
    Cerr("***** castor-datafileConversionEx :: A problem occurred while logging command line arguments !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Output progression options
  cout << std::fixed;
  cout << std::setprecision(2);

  // ============================================================================================================
  // ScannerManager object initialization
  // (Basic initialization to locate scanner geometry files)
  // ============================================================================================================

  sScannerManager *p_scannermanager = sScannerManager::GetInstance();
  p_scannermanager->SetVerbose(vb);

  // Check if the scanner exists and get the name from ScannerManager
  scanner_alias = (scanner_alias.find(OS_SEP)) ? scanner_alias.substr(scanner_alias.find_last_of(OS_SEP) + 1) : scanner_alias;

  if (p_scannermanager->FindScannerSystem(scanner_alias))
  {
    Cerr("**** castor-datafileConversionEx :: A problem occurred while searching for scanner system !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Initialize the scanner system object
  p_scannermanager->BuildScannerObject();
  p_scannermanager->GetScannerObject()->Instantiate(true);
  p_scannermanager->CheckParameters();
  p_scannermanager->Initialize();
  cout << endl;
  cout << "----------------------------------------" << endl;
  cout << "Description of the scanner system object" << endl;
  cout << "----------------------------------------" << endl;
  p_scannermanager->Describe();
  cout << endl;
  
  // Build output file names
  path_to_data_file = path_to_out_file + ".cdf";
  path_to_header_file = path_to_out_file + ".cdh";

  // ============================================================================================================
  // Instanciate/Initialize CASToR DataFile
  // (Enable/Disable optionnal field to be written)
  // ============================================================================================================

  // Instantiate & Initialize oImageDimensionsAndQuantification object, required for datafile generation (number of threads)
  oImageDimensionsAndQuantification *p_ID = new oImageDimensionsAndQuantification();
  p_ID->SetDefault();
  p_ID->SetVerbose(vb);
  // p_ID->CheckParameters();
  // p_ID->SetNbVoxX();
  p_ID->SetNbBeds(1);
  p_ID->Initialize();
  // p_ID->SetNbThreads("4");

  // Instanciate & Initialize iDataFilePET and Event objects
  iDataFilePET *Out_data_file = NULL;
  vEvent *Event = NULL;

  // ===> This variable set the maximum number of LORs contained in one event
  // ===> in the whole dataset (for compression purposes)
  // ===> This depends on the dataset to convert
  uint16_t max_nb_lines_per_event = 1;

  // Create datafile and set variables
  Out_data_file = new iDataFilePET();

  // Set the image oImageDimensionsAndQuantification object (required for initialization)
  Out_data_file->SetImageDimensionsAndQuantification(p_ID);

  // Set path to header datafile
  Out_data_file->SetHeaderDataFileName(path_to_header_file);

  // Set verbosity
  Out_data_file->SetVerbose(vb);

  // Set data mode (list/hitogram)
  Out_data_file->SetDataMode(data_mode);

  // Set modality
  Out_data_file->SetDataType(TYPE_PET);

  // Set maximum number of LORs contained in one event
  Out_data_file->SetMaxNumberOfLinesPerEvent(max_nb_lines_per_event);

  // Set calibration factor
  Out_data_file->SetCalibrationFactor(calib_factor);
  
  // Set isotope
  if (!istp_alias.empty())
  {
    if (p_ID->SetPETIsotope(0, istp_alias))
    {
      Cerr("**** castor-datafileConversionEx :: A problem occurred while checking isotope name !" << endl);
      Exit(EXIT_FAILURE);
    }
    Out_data_file->SetIsotope(istp_alias);
  }

  // Allocate the Event structure (containing all information for a listmode
  // or histogram event, depending on data mode)
  if (data_mode == MODE_HISTOGRAM)
  {
    // Instanciate histogram Event
    Event = new iEventHistoPET();
    Event->SetNbLines(max_nb_lines_per_event);
    Event->AllocateID();
  }
  else
  {
    // Instanciate list-mode Event
    Event = new iEventListPET();
    Event->SetNbLines(max_nb_lines_per_event);
    Event->AllocateID();
  }

  // ============================================================================================================
  // ===> Read sinogram correction file(s) (if any) and set datafile flag for each correction
  // ============================================================================================================

  // ===> Compute sinogram dimensions.
  // ===> Values must be adjusted to the dataset to convert.
  uint32_t nbins = 1;
  uint32_t nangles = 1;
  uint32_t nsinos = 1;
  uint32_t nb_cor_factors = nbins * nangles * nsinos;

  // --- Normalization correction factors ---

  // FLTNBDATA is either a float (default) or double type.
  // Defined in /management/gVariables.hh
  FLTNBDATA *p_norm = NULL;
  bool norm_flag = false;

  if (!path_to_norm_file.empty())
  {
    if (vb >= 2)
      Cout("--> Reading the norm correction factors file = " << path_to_norm_file << "..." << endl);
    ifstream ifile(path_to_norm_file.c_str(), ios::binary | ios::in);

    if (!ifile)
    {
      Cerr("***** castor-datafileConversionEx ::  Error reading normalization factor file: " << path_to_norm_file << " !" << endl);
      Exit(EXIT_FAILURE);
    }
    else
    {
      // Allocate buffer storing the normalisation
      p_norm = new FLTNBDATA[nb_cor_factors * sizeof(FLTNBDATA)];
      ifile.read(reinterpret_cast<char *>(p_norm), nb_cor_factors * sizeof(FLTNBDATA));
      ifile.close();
    }

    // ===> Set normalization correction flag in datafile objet
    // ===> (required to indicate that normalization correction will be provided)
    Out_data_file->SetNormCorrectionFlagOn();
    norm_flag = true;
  }

  // --- Scatter correction factors file ---
  FLTNBDATA *p_scat = NULL;
  bool scat_flag = false;

  if (!path_to_scat_file.empty())
  {
    if (vb >= 2)
      Cout("--> Reading the scatter correction factors file = " << path_to_scat_file << "..." << endl);
    ifstream ifile(path_to_scat_file.c_str(), ios::binary | ios::in);

    if (!ifile)
    {
      Cerr("***** castor-datafileConversionEx ::  Error reading scatter factor file: " << path_to_scat_file << " !" << endl);
      Exit(EXIT_FAILURE);
    }
    else
    {
      // Allocating buffer storing the scatter
      p_scat = new FLTNBDATA[nb_cor_factors * sizeof(FLTNBDATA)];
      ifile.read(reinterpret_cast<char *>(p_scat), nb_cor_factors * sizeof(FLTNBDATA));
      ifile.close();
    }

    // ===> Set scatter correction flag in datafile objet
    // ===> (required to indicate that scatter correction will be provided)
    Out_data_file->SetScatterCorrectionFlagOn();
    scat_flag = true;
  }

  // --- Random correction factors file ---
  FLTNBDATA *p_rdm = NULL;
  bool rdm_flag = false;

  if (!path_to_rdm_file.empty())
  {
    if (vb >= 2)
      Cout("--> Reading the random correction factors file = " << path_to_rdm_file << "..." << endl);
    ifstream ifile(path_to_rdm_file.c_str(), ios::binary | ios::in);

    if (!ifile)
    {
      Cerr("***** castor-datafileConversionEx ::  Error reading random correction factor file: " << path_to_rdm_file << " !" << endl);
      Exit(EXIT_FAILURE);
    }
    else
    {
      // Allocating buffer storing the random
      p_rdm = new FLTNBDATA[nb_cor_factors * sizeof(FLTNBDATA)];
      ifile.read(reinterpret_cast<char *>(p_rdm), nb_cor_factors * sizeof(FLTNBDATA));
      ifile.close();
    }

    // ===> Set random correction flag in datafile objet
    // ===> (required to indicate that random correction will be provided)
    Out_data_file->SetRandomCorrectionFlagOn();
    rdm_flag = true;
  }

  // --- Attenuation correction factors file ---
  FLTNBDATA *p_acf = NULL;
  bool acf_flag = false;

  if (!path_to_acf_file.empty())
  {
    if (vb >= 2)
      Cout("--> Reading the attenuation correction factors file = " << path_to_acf_file << "..." << endl);
    ifstream ifile(path_to_acf_file.c_str(), ios::binary | ios::in);

    if (!ifile)
    {
      Cerr("***** castor-datafileConversionEx ::  Error reading attenuation correction factor file: " << path_to_acf_file << " !" << endl);
      Exit(EXIT_FAILURE);
    }
    else
    {
      // Allocating buffer storing the acf
      p_acf = new FLTNBDATA[nb_cor_factors * sizeof(FLTNBDATA)];
      ifile.read(reinterpret_cast<char *>(p_acf), nb_cor_factors * sizeof(FLTNBDATA));
      ifile.close();
    }

    // ===> Set attenuation correction flag in datafile objet
    // ===> (required to indicate that attenuation correction will be provided)
    Out_data_file->SetAtnCorrectionFlagOn();
    acf_flag = true;
  }

  // Initialize datafile object
  Out_data_file->PROJ_InitFile();

  // Compute the event size relatively to the mandatory
  // and optionnal (correction factors) elements
  Out_data_file->ComputeSizeEvent();

  // Init DataFile
  Out_data_file->PrepareDataFile();


  // ============================================================================================================
  // *********************************************** CONVERSION *************************************************
  // ============================================================================================================
  Cout(" --- Start conversion of datafile : " << path_to_input_file << endl
                                              << "     CASToR output header datafile: " << path_to_header_file << endl
                                              << "     CASToR output binary datafile: " << path_to_data_file << endl
                                              << endl);

  // ===> Open the input datafile
  ifstream infile_input(path_to_input_file);
	string line;
 
  // ===> Check if the input datafile is open
  if (!infile_input.is_open())
  {
    Cerr("***** castor-datafileConversionEx ::  Error opening input datafile: " << path_to_input_file << " !" << endl);
    Exit(EXIT_FAILURE);
  }


  // ============================================================================================================
  // Init variables
  // ============================================================================================================

  // Index for progression printing
  uint32_t printing_index = 1;

  // Skip the header line
	getline(infile_input, line); 

  // ===> Recover the total number of events (number of coincidences or histogram bins)
  uint32_t nbEvents = 0;
  // Save the values of times and ids for each event
  vector<uint32_t> times;
  vector<uint32_t> ids1;
  vector<uint32_t> ids2;
  while (std::getline(infile_input, line)){
    ++nbEvents;
    istringstream iss(line);
    float time;
    uint32_t id1, id2;
    // // If separated by comma
    // char comma; // comma character
    // if (!(iss >> time >> comma >> id1 >> comma >> id2)) { break; } // error
    if (!(iss >> time >> id1 >> id2)) { break; } // error
    times.push_back(time);
    ids1.push_back(id1);
    ids2.push_back(id2);
  }
  infile_input.close();
  cout << "nbEvents = " << nbEvents << endl << endl;

  // ===> Acquisition duration and start time in seconds
  FLTNB start_time_sec = times[0] / 1000;
  FLTNB stop_time_sec = times.back() / 1000;
  FLTNB duration_sec = stop_time_sec - start_time_sec;

  // ===> Variables to recover indices of the two crystals for each potential
  // ===> LORs contributing to the event
  uint32_t *castor_id1 = new uint32_t[max_nb_lines_per_event];
  uint32_t *castor_id2 = new uint32_t[max_nb_lines_per_event];


  // Just for progression output
  uint32_t printing_ratio = (nbEvents > 1000) ? 1000 : nbEvents / 10;

  // Loop on the data elements (histogram bins or LORs)
  for (uint32_t i = 0; i < nbEvents; i++)
  {
    // Define IDs variables
    uint32_t crystalID1 = ids1[i], crystalID2 = ids2[i];

    // ===> Recover event value (list-mode, or dynamic histogram dataset)
    uint32_t time_ms = times[i];
    stop_time_sec = (FLTNB)time_ms / 1000;

    // ===> Recover the number of lines in the event (if compression).
    int nb_lines_in_event = 1; // No compression

    // Loop on the LORs contained in the event
    for (int l = 0; l < nb_lines_in_event; l++)
    {
      // Recover the CASToR ID position of the two crystals according to
      // the LUT corresponding to this scanner
      castor_id1[l] = crystalID1;
      castor_id2[l] = crystalID2;
    }
    // -----------------------------------------------------------------
    // -----------------------------------------------------------------

    // --- Histo-mode event ---
    if (is_histogram_flag) // histogram mode event
    {
      // Set the number of lines in this event
      Event->SetNbLines(nb_lines_in_event);

      // Set CASToR IDs in the event structure for each contributing line
      for (int l = 0; l < nb_lines_in_event; l++)
      {
        Event->SetID1(l, castor_id1[l]); // 1st argument : line, 2nd argument : index
        Event->SetID2(l, castor_id2[l]); // 1st argument : line, 2nd argument : index
      }

      // Set information related to each TOF bin
      int nb_tof_bins = 1;
      for (int b = 0; b < nb_tof_bins; b++)
      {
        // ===> Recover and set event value (histogram dataset)
        FLTNB event_value = 1.;
        Event->SetEventValue(b, event_value); // 1st argument : TOF bin, 2nd argument : event value

        // ===> Set scatter correction factor (if any)
        // ===> Require computation of the scatter correction value for
        // ===> this LOR from the factors loaded in the "p_scat" buffer
        if (scat_flag)
        {
          FLTNB scat_correction_coeff = 0.; // computed from p_scat
          ((iEventHistoPET *)Event)->SetScatterRate(b, scat_correction_coeff);
        }
      }

      // ===> Set timestamp of the event in ms (global timestamp for histogram)
      // ===> i.e : 0, or timestamp of a dynamic frame
      Event->SetTimeInMs(time_ms); // uint32_t

      // ===> Set random correction factor (if any)
      // ===> Require computation of the random correction value for
      // ===> this LOR from the factors loaded in the "p_rdm" buffer
      if (rdm_flag)
      {
        FLTNB rdm_correction_coeff = 0.; // computed from p_rdm
        ((iEventPET *)Event)->SetRandomRate(rdm_correction_coeff);
      }

      // ===> Set normalization correction factor (if any)
      // ===> Require computation of the norm correction value for
      // ===> this LOR from the factors loaded in the "p_norm" buffer
      if (norm_flag)
      {
        FLTNB norm_correction_coeff = 0.; // computed from p_norm
        ((iEventPET *)Event)->SetNormalizationFactor(norm_correction_coeff);
      }

      // ===> Set attenuation correction factor (if any)
      // ===> Require computation of the atn correction value for
      // ===> this LOR from the factors loaded in the "p_acf" buffer
      if (acf_flag)
      {
        FLTNB atn_correction_coeff = 0.; // computed from p_acf
        ((iEventPET *)Event)->SetAttenuationCorrectionFactor(atn_correction_coeff);
      }

      // Write event in the datafile
      // ('0' argument is just thread-related. No change required)
      Out_data_file->WriteEvent(Event, 0);
    }

    // -----------------------------------------------------------------
    // --- List-mode event ---
    else
    {
      // Set the number of lines in this event
      Event->SetNbLines(nb_lines_in_event);

      // Set CASToR ID in the event structure
      for (int l = 0; l < nb_lines_in_event; l++)
      {
        Event->SetID1(l, castor_id1[l]); // 1st argument : line, 2nd argument :index
        Event->SetID2(l, castor_id2[l]); // 1st argument : line, 2nd argument :index
      }

      // ===> Set scatter correction factor (if any)
      // ===> Require computation of the scatter correction value for
      // ===> this LOR from the factors loaded in the "p_scat" buffer
      if (scat_flag)
      {
        FLTNB scat_correction_coeff = 0.; // computed from p_scat
        // First argument '0' by default for list-mode (it is used for histograms with TOF bins).
        ((iEventListPET *)Event)->SetScatterRate(0, scat_correction_coeff);
      }

      // ===> Set timestamp of the event in ms
      Event->SetTimeInMs(time_ms); // uint32_t

      // ===> Set random correction factor (if any)
      // ===> Require computation of the random correction value for
      // ===> this LOR from the factors loaded in the "p_rdm" buffer
      if (rdm_flag)
      {
        FLTNB rdm_correction_coeff = 0.; // computed from p_rdm
        ((iEventPET *)Event)->SetRandomRate(rdm_correction_coeff);
      }

      // ===> Set normalization correction factor (if any)
      // ===> Require computation of the norm correction value for
      // ===> this LOR from the factors loaded in the "p_norm" buffer
      if (norm_flag)
      {
        FLTNB norm_correction_coeff = 0.; // computed from p_norm
        ((iEventPET *)Event)->SetNormalizationFactor(norm_correction_coeff);
      }

      // ===> Set attenuation correction factor (if any)
      // ===> Require computation of the atn correction value for
      // ===> this LOR from the factors loaded in the "p_acf" buffer
      if (acf_flag)
      {
        FLTNB atn_correction_coeff = 0.; // computed from p_acf
        ((iEventPET *)Event)->SetAttenuationCorrectionFactor(atn_correction_coeff);
      }

      // Write event in the datafile
      // ('0' argument is just thread-related. No change required)
      Out_data_file->WriteEvent(Event, 0);
    }

    // Progression output
    if (printing_index % printing_ratio == 0)
    {
      FLTNB percent = (((FLTNB)(printing_index + 1)) / ((FLTNB)nbEvents)) * ((FLTNB)100);
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
           << percent << " %                    ";
    }
    printing_index++;
  }

  Cout("The simulated dataset contained " << nbEvents << " coincidences / bins" << endl);

  // Set duration and nb of events in datafile
  Out_data_file->SetStartTime(start_time_sec);
  Out_data_file->SetDuration(duration_sec);
  Out_data_file->SetNbEvents(nbEvents);

  // Write datafile header
  Out_data_file->WriteHeader();

  // Write raw datafile
  Out_data_file->PROJ_WriteData();

  // Delete tmp writing file (if any)
  Out_data_file->PROJ_DeleteTmpDataFile();

  cout << endl;

// ============================================================================================================
// End
// ============================================================================================================

// Free Root objects
#ifdef CASTOR_ROOT
  delete Coincidences;
  delete Tfile_root;
  delete Tapp;
#endif

  delete[] castor_id1;
  delete[] castor_id2;

  if (p_acf != NULL)
    delete[] p_acf;
  if (p_rdm != NULL)
    delete[] p_rdm;
  if (p_scat != NULL)
    delete[] p_scat;
  if (p_norm != NULL)
    delete[] p_norm;

  // Delete objects
  if (Event)
  {
    if (is_histogram_flag)
      delete (iEventHistoPET *)Event;
    else
      delete (iEventListPET *)Event;
  }

  delete Out_data_file;
  delete p_ID;

  return 0;
}
