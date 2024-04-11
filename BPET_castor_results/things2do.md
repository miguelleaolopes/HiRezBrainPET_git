# Progress of RPC reconstruction

## Scanner Geometry

Have a system, with crystal with a size of 1mm transaxial and axial.
After done, redo for 1/4 of size?? But this will make x4x4 more elements, and its already pretty big (286k elets).

## Modify python script for the compiling

create a batch script dependent on the system?

## Do a ReadMe Documentation

Create a ReadMe file with the instructions to compile and run the code, as well as the parameters to change and the options to see.

- Start with the Intersect_LOR.py to transform the datafile given into fractions with the time and coord
- Use the python/combine_files.py to combine the fraction data into a single txt file LOR_all_points_time.txt
- Check if the scanner wanted is in the config/scanner folder
  - If not, create with the castor-PETScannerLutEx.exe
  - There is a Release/lut_explorer.bat to help confirm the scanner
- Use the Release/txt_id_conversion.bat to convert the LOR coord to the id of the crystal in the scanner geometry
- Use the Release/data_conversion.bat to convert the txt file into CASToR format
- Given the final data, one can start using the main program of reconstruction, castor_recon.exe
  - The parameters are all combine in a batch file to run, run{...}_win{}.bat

Prerequisites - of the list that i have, e.g. python, cmake, mpi, etc

Instalation - of the requirements.txt

Usage - of the scripts, with the parameters and options. And of the python file. 2 different sections.

## Other stuff

python computeDerenzoValleyToPeak.py -f example/1test_it6.hdr -c example/config.json -m mean --showTriangPos --showSpotsPos --showLinesProfile --imageShown 2 --showVprHistos --saveResults

