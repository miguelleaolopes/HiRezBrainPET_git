/*!
    \file
    \ingroup utils
    \brief This program is used to convert a datafile txt with the positions of crystal to IDs of crystals, taken into account the LUT of the detector.
    \details The program reads a CASToR Look-UP-Table (LUT) file and a datafile txt with the positions of crystal
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "iDataFilePET.hh"
// #include "iDataFileSPECT.hh"
// #include "iDataFileCT.hh"
#include "sScannerManager.hh"
#include "iScannerPET.hh"

/**
 * @defgroup INPUT_FILE_TYPE input file type
 *
 *    \brief Input file type for castor-scannerLUTExplorer (lut or geom) \n
 *           Defined in castor-scannerLUTExplorer.cc
 * @{
 */
/** Constant corresponding to an undefined type (default) (=-1) */
#define IF_TYPE_UNKNOWN -1
/** Constant corresponding to a Look-Up-Table file (.hscan/.lut) (=0) */
#define IF_TYPE_LUT 0
/** Constant corresponding to a generic geometry file (.geom )file (=1) */
#define IF_TYPE_GEOM 1
/** @} */

// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================
//                                                        H E L P     F U N C T I O N S
// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================

/*!
  \fn      ShowHelp()
  \brief   Display main command line options for castor-recon
*/
void ShowHelp()
{
  // Show help
  cout << endl;
  cout << "Usage: castor-txtConversionCrystalsID -txt datadfile.txt -sf config/scanner/scanfile.hscan -o txtoutput.txt" << endl;
  cout << endl;
  cout << "This program can be used to convert a datafile txt with the positions of crystal to IDs of crystals, taken into account the LUT of the detector." << endl;
  cout << endl;
  cout << endl;
  cout << "[Main Settings]:" << endl;
  cout << "  -txt <filename>      : datafile txt with the positions of crystal" << endl;
  cout << "  -sf <filename>       : give the path to the header of lut scanfile" << endl;
  cout << "  -o <filename>        : give the path to the output datafile txt with the IDs of crystals" << endl;
  cout << "                         (in relation to the txt location, if just name, saves in same as txt path)" << endl;
  cout << endl;
  cout << "[Miscellaneous Options]:" << endl;
  cout << "  -h                   : Print out this help page" << endl;
  cout << "  -vb value            : Verbose mode" << endl;
  cout << endl;
}

// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================
//                                                        M A I N     P R O G R A M
// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================

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
  {
    ShowHelp();
    Exit(EXIT_SUCCESS);
  }

  // =============================================================================================================================================
  // Define a few parameters
  // =============================================================================================================================================

  // String to store the path of the input filename
  string path_to_txt_filename = "";

  // String to store the path of the scanner filename
  string scanner_alias = "";

  // String to store the path of the output filename
  string path_to_output_filename = "";

  // General verbose level
  int verbose = 0;
  // Fald for input file type (0 for LUT, 1 for GEOM)
  int input_file_type = IF_TYPE_LUT;

  for (int i = 1; i < argc; i++)
  {
    // Get the option as a string
    string option = (string)argv[i];
    // Show help
    if (option == "-h")
    {
      ShowHelp();
      Exit(EXIT_SUCCESS);
    }
    // General verbose level
    else if (option == "-vb")
    {
      if (i >= argc - 1)
      {
        Cerr("***** castor-txtConversionCrystalsID() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i + 1], &verbose))
      {
        Cerr("***** castor-txtConversionCrystalsID() -> Exception when trying to read provided verbosity level '" << verbose << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      verbose = atoi(argv[i + 1]);
      i++;
    }
    // ScannerFile
    else if (option == "-sf")
    {
      if (i >= argc - 1)
      {
        Cerr("***** castor-txtConversionCrystalsID() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      scanner_alias = (string)argv[i + 1];
      i++;
    }
    // DataFile
    else if (option == "-txt")
    {
      if (i >= argc - 1)
      {
        Cerr("***** castor-txtConversionCrystalsID() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_txt_filename = (string)argv[i + 1];
      i++;
    }
    // OutputFile
    else if (option == "-o")
    {
      if (i >= argc - 1)
      {
        Cerr("***** castor-txtConversionCrystalsID() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_output_filename = (string)argv[i + 1];
      i++;
    }
    // Unknown option
    else
    {
      Cerr("***** castor-txtConversionCrystalsID() -> Unknown option: " << option << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // =============================================================================================================================================
  // Check that all mandatory options have been provided and are valid
  // =============================================================================================================================================
  if (path_to_txt_filename == "")
  {
    Cerr("***** castor-txtConversionCrystalsID() -> No datafile provided" << endl);
    Exit(EXIT_FAILURE);
  }
  if (scanner_alias == "")
  {
    Cerr("***** castor-txtConversionCrystalsID() -> No scanfile provided" << endl);
    Exit(EXIT_FAILURE);
  }
  if (path_to_output_filename == "")
  {
    Cerr("***** castor-txtConversionCrystalsID() -> No outputfile provided" << endl);
    Exit(EXIT_FAILURE);
  }
  if (verbose < 0)
  {
    Cerr("***** castor-txtConversionCrystalsID() -> Verbose less than 1 has no sense for a program used to solely print out information on screen !" << endl);
    Exit(EXIT_FAILURE);
  }

  // =============================================================================================================================================
  // Initializations
  // =============================================================================================================================================

  // Verbose (we know it is at least 1 so we don't check it)
  Cout("==============================================================" << endl);
  Cout("castor-txtConversionCrystalsID() -> Initialization starts" << endl);

  // Get user endianness (interfile I/O)
  GetUserEndianness();

  // =============================================================================================================================================
  // Create sScannerManager and Read the LUT file
  // =============================================================================================================================================

  sScannerManager *p_scannerManager = sScannerManager::GetInstance();
  p_scannerManager->SetVerbose(verbose);

  // Check if the scanner exists and get the name from ScannerManager
  scanner_alias = (scanner_alias.find(OS_SEP)) ? scanner_alias.substr(scanner_alias.find_last_of(OS_SEP) + 1) : scanner_alias;
  Cout("Scanner alias: " << scanner_alias << endl);

  // Read the LUT file
  if (p_scannerManager->FindScannerSystem(scanner_alias))
  {
    Cerr("**** castor-datafileConversionEx :: A problem occurred while searching for scanner system !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Initialize scanner object
  p_scannerManager->BuildScannerObject();
  p_scannerManager->InstantiateScanner();
  p_scannerManager->CheckParameters();
  p_scannerManager->Initialize();
  p_scannerManager->GetScannerObject()->LoadLUT();
  cout << endl;
  cout << "----------------------------------------" << endl;
  cout << "Description of the scanner system object" << endl;
  cout << "----------------------------------------" << endl;
  p_scannerManager->Describe();
  cout << endl;

  // =============================================================================================================================================
  // Get number of elements in the system
  int nb_elts = p_scannerManager->GetSystemNbElts();

  // Get number of rings axially
  int nb_rings_lyr;
  ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of rings in scanner", &nb_rings_lyr, 1, KEYWORD_OPTIONAL);

  // Get number of rsectors - this assumes that there are no modules and submodules (mod=submod=1)
  int nb_rsectors = nb_elts / nb_rings_lyr;

  // -------------------------------------------------------------------------------------------------------
  // Create list of crystals positions and orientations 1 and 2
  FLTNB **crystal_positions1 = new FLTNB *[nb_elts];
  FLTNB **crystal_positions2 = new FLTNB *[nb_elts];
  FLTNB **crystal_orientations1 = new FLTNB *[nb_elts];
  FLTNB **crystal_orientations2 = new FLTNB *[nb_elts];

  for (int i = 0; i < nb_elts; i++)
  {
    crystal_positions1[i] = new FLTNB[3];
    crystal_positions2[i] = new FLTNB[3];
    crystal_orientations1[i] = new FLTNB[3];
    crystal_orientations2[i] = new FLTNB[3];
  }

  // Get index start and stop
  int64_t index_start = 0;
  int64_t index_stop = p_scannerManager->GetSystemNbElts();

  // Loop through each element and get its position and orientation
  for (int index = index_start; index < index_stop; index++)
  {
    // Call the function to get positions and orientations for the current index
    int result = p_scannerManager->GetScannerObject()->GetPositionsAndOrientations(0, index, crystal_positions1[index], crystal_positions2[index], crystal_orientations1[index], crystal_orientations2[index]);

    if (result != 0)
    {
      Cerr("**** castor-datafileConversionEx :: A problem occurred while getting positions and orientations for index " << index << " !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  cout << "nb_elts: " << nb_elts << endl;
  cout << "nb_rings_lyr: " << nb_rings_lyr << endl;
  cout << "nb_rsectors: " << nb_rsectors << endl;
  cout << "------------------------------------------------------" << endl;

  // =============================================================================================================================================
  // Read txt data file (path_to_txt_filename), convert the IDs, and write the output to a new txt file (path_to_output_filename)
  // =============================================================================================================================================

  // Define output filename 
  // Get input location of the txt file
  string path_to_txt_filename_no_ext = path_to_txt_filename.substr(0, path_to_txt_filename.find_last_of(OS_SEP));
  // Add the output filename
  path_to_output_filename = path_to_txt_filename_no_ext + OS_SEP + path_to_output_filename;

  // Open the output file
  ofstream output_file;
  output_file.open(path_to_output_filename);
  // Check if the file is open
  if(!output_file.is_open())
  {
    Cerr("**** castor-datafileConversionEx :: A problem occurred while opening the output file !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Write the header to the file
  output_file << "# time (msec) id1 id2" << std::endl;

  // Variables to read the txt file
  ifstream infile(path_to_txt_filename);
  string line;

  int totalLines = 0;
  // Count the total number of lines in the file
  // Skip first line
  Cout("Counting the total number of lines in the file..." << endl);
  getline(infile, line);
  while (getline(infile, line)) {
      totalLines++;
  }
  Cout("Total number of lines: " << totalLines << endl);
  Cout("==============================================================" << endl);

  // Reset the file stream and progress counters
  infile.clear();
  infile.seekg(0);

  // Skip first line
  getline(infile, line);
  // For testing only read first 3 lines
  int line_count = 0;
  
  // Index for progression printing
  uint32_t printing_index = 0;
  // Just for progression output
  // uint32_t printing_ratio = (totalLines > 1000) ? 1000 : totalLines / 10;

  // -----------------------------------------------------------------
  // Loop through each line
  while (getline(infile, line))
  {
    // Create a stringstream of the current line and parse the comma-separated values
    istringstream iss(line);
    FLTNB time, x1, y1, z1, x2, y2, z2;

    char comma; // To handle the comma separator
    iss >> time >> comma >> x1 >> comma >> y1 >> comma >> z1 >> comma >> x2 >> comma >> y2 >> comma >> z2;

    // Check if parsing was successful
    if (iss.fail())
    {
      cout << "**** castor-datafileConversionEx :: A problem occurred while reading the txt file !" << endl;
      Exit(EXIT_FAILURE);
    }

    // Finding nearest crystal in the scanner
    int min_index1 = -1;
    int min_index1_z = -1;
    int min_index2 = -1;
    int min_index2_z = -1;

    #undef max
    FLTNB min_distance1 = std::numeric_limits<FLTNB>::max();
    FLTNB min_distance2 = std::numeric_limits<FLTNB>::max();

    // Loop through each crystal in the scanner, deprecated, its a slower process, bigger cycle
    // for(int i = 0; i < nb_elts; i++) {
    //   FLTNB distance1 = sqrt(pow(x1 - crystal_positions2[i][0], 2) + pow(y1 - crystal_positions2[i][1], 2) + pow(z1 - crystal_positions2[i][2], 2));
    //   FLTNB distance2 = sqrt(pow(x2 - crystal_positions2[i][0], 2) + pow(y2 - crystal_positions2[i][1], 2) + pow(z2 - crystal_positions2[i][2], 2));
    //   if (min_distance1 > distance1) {
    //     min_distance1 = distance1;
    //     min_index1 = i;
    //   }
    //   if (min_distance2 > distance2) {
    //     min_distance2 = distance2;
    //     min_index2 = i;
    //   }
    // }

    // Other way to do it - considering the cylindrical symmetry of the scanner

    // Look first for the nearest crystal axially in the z direction
    for (int i = 0; i < nb_elts; i += nb_rsectors)
    {
      FLTNB distance1 = abs(z1 - crystal_positions2[i][2]);
      FLTNB distance2 = abs(z2 - crystal_positions2[i][2]);
      if (min_distance1 > distance1)
      {
        min_distance1 = distance1;
        min_index1_z = i;
      }
      if (min_distance2 > distance2)
      {
        min_distance2 = distance2;
        min_index2_z = i;
      }
    }

    min_distance1 = std::numeric_limits<FLTNB>::max();
    min_distance2 = std::numeric_limits<FLTNB>::max();
    // Look for the nearest crystal in the transverse plane
    for (int i = 0; i < nb_rsectors; i++)
    {
      FLTNB distance1 = sqrt(pow(x1 - crystal_positions2[min_index1_z + i][0], 2) + pow(y1 - crystal_positions2[min_index1_z + i][1], 2));
      FLTNB distance2 = sqrt(pow(x2 - crystal_positions2[min_index2_z + i][0], 2) + pow(y2 - crystal_positions2[min_index2_z + i][1], 2));
      if (min_distance1 > distance1)
      {
        min_distance1 = distance1;
        min_index1 = min_index1_z + i;
      }
      if (min_distance2 > distance2)
      {
        min_distance2 = distance2;
        min_index2 = min_index2_z + i;
      }
    }

    // Save the results in the output file
    output_file << time << " " << min_index1 << " " << min_index2 << std::endl;

    // Prints the results to check
    // cout << "Time: " << time << endl;
    // cout << "Position 1: " << x1 << " " << y1 << " " << z1 << endl;
    // cout << "Nearest crystal index1: " << min_index1 << endl;
    // cout << "Nearest crystal position 1: " << crystal_positions2[min_index1][0] << " " << crystal_positions2[min_index1][1] << " " << crystal_positions2[min_index1][2] << endl;
    // cout << "Position 2: " << x2 << " " << y2 << " " << z2 << endl;
    // cout << "Nearest crystal index2: " << min_index2 << endl;
    // cout << "Nearest crystal position 2: " << crystal_positions2[min_index2][0] << " " << crystal_positions2[min_index2][1] << " " << crystal_positions2[min_index2][2] << endl;
    // cout << "------------------------------------------------------" << endl;
    // cout << "Min distance 1: " << min_distance1 << endl;
    // cout << "Min distance 2: " << min_distance2 << endl;
    // cout << "------------------------------------------------------" << endl;

    // At the end of the loop just to test a few lines:
    // line_count += 1;
    // if (line_count == 5)
    // {
    //   break;
    // }
    // Progression output
    // if (printing_index % printing_ratio == 0)
    // {
    //   FLTNB percent = (((FLTNB)(printing_index + 1)) / ((FLTNB)totalLines)) * ((FLTNB)100);
    //   cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
    //        << percent << " %                    ";
    // }
    
    // Update the progress percentage
    float progress = static_cast<float>(printing_index) / totalLines * 100;
    
    std::cout << "\rProgress: " << std::fixed << std::setprecision(2) << progress << " %";
    
    std::cout.flush();  // Flush the output buffer to update the line immediately
    printing_index++;
  }

  // End of the loop
  Cout(endl); 
  Cout("==============================================================" << endl);
  Cout(":::::::::::::: Conversion done ! ::::::::::::::" << endl);

  // Close the output file and the input file
  output_file.close();
  infile.close();

  // Deallocate the memory
  for (int i = 0; i < nb_elts; i++)
  {
    delete[] crystal_positions1[i];
    delete[] crystal_positions2[i];
    delete[] crystal_orientations1[i];
    delete[] crystal_orientations2[i];
  }
  delete[] crystal_positions1;
  delete[] crystal_positions2;
  delete[] crystal_orientations1;
  delete[] crystal_orientations2;

  // ============================================================================================================
  // End
  // ============================================================================================================

  return EXIT_SUCCESS;
}