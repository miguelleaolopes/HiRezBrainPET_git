@echo off

@REM This script converts the data from the RPC_HIREZBRAINPET scanner to the CASToR format.
@REM The data is stored in the file intersect_points_with_time.txt

@REM Define variables for the command line
set command=.\castor-datafileConversionEx.exe
set isotope=-ist F-18
set verbose=-vb 2
@REM set calibration_fac = -cf 1.0

set scanner=-s ..\..\config\scanner\HiRezBrainPET_rsmall_Geom 
set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_orig
set datafile=-il ..\Files\LOR_all_time_ID_Geom_rsmall_orig.txt

@REM set datafile=-il ..\Files\LOR_all_time_ID_z2s.txt
@REM set output=-o ..\Files\CASToR_data\RPC_all_data_z2s
@REM set scanner=-s ..\..\config\scanner\RPC_HiRezBrainPET_Geom_z2s

@REM set datafile=-il ..\Files\LOR_all_time_ID_small.txt
@REM set output=-o ..\Files\CASToR_data\RPC_all_data_small
@REM set scanner=-s ..\..\config\scanner\RPC_HiRezBrainPET_Geom_small

@REM set datafile=-il ..\Files\LOR_all_time_ID_big.txt
@REM set output=-o ..\Files\CASToR_data\RPC_all_data_big
@REM set scanner=-s ..\..\config\scanner\RPC_HiRezBrainPET_Geom_big

echo ==================================================================================
echo Conversion of RPC_HIREZBRAINPET data.
echo ==================================================================================
%command% %datafile% %output% %scanner% %verbose% %isotope% 
