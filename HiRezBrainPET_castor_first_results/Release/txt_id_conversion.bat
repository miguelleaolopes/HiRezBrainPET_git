@echo off

@REM This script converts the txt data with the positions of the LOR (t,x1,y1,z1,x2,y2,z2) to crystal ID.
@REM The data is stored in the file LOR_all_points_time.txt

@REM Define variables for the command line
set command=.\castor-txtConversionCrystalsID.exe
set verbose=-vb 3
set txtfile=-txt ..\Files\LOR_all_points_time.txt

set output=-o LOR_all_time_ID_Geom_rsmall_orig.txt 
set scanner=-sf ..\config\scanner\HiRezBrainPET_rsmall_Geom

@REM set output=-o LOR_all_time_ID_z2s.txt 
@REM set scanner=-sf ..\config\scanner\RPC_HiRezBrainPET_Geom_z2s

@REM set output=-o LOR_all_time_ID_small.txt
@REM set scanner=-sf ..\config\scanner\RPC_HiRezBrainPET_Geom_small

@REM set output=-o LOR_all_time_ID_big.txt
@REM set scanner=-sf ..\config\scanner\RPC_HiRezBrainPET_Geom_big


echo ==================================================================================
echo Conversion of TXT file with crystal ID.
echo ==================================================================================
%command% %txtfile% %output% %scanner% %verbose%
