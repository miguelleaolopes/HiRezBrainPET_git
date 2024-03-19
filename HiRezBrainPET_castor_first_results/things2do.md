# Progress of RPC reconstruction

## Scanner Geometry

Have a system, with crystal with a size of 1mm transaxial and axial.
After done, redo for 1/4 of size?? But this will make x4x4 more elements, and its already pretty big (286k elets).

## New bat file with only the sets

Create a new bat file with the set commands for the new pc, test if the exe already works in the new pc without compiling again.

## Create a python script for the compiling

Create a python script that makes a gui to chose the parameters of the compiling, and then compile the code by creating a batch script (system versions?).

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

## Add the option to add the scanner either in the config/scanner folder or a specific path with the InitScannerWithFile function of sScannerManager.cc

## Change CmakeLists.txt to include the modified names for the files for the HirRezBrainPET
