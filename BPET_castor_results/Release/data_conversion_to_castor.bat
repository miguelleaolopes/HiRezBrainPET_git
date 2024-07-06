@echo off

@REM This script converts the data from the RPC_HIREZBRAINPET scanner to the CASToR format.
@REM The data is stored in the file intersect_points_with_time.txt

@REM Define variables for the command line
set command=.\castor-datafileConversionEx_BPET.exe
set isotope=-ist F-18
set verbose=-vb 2
@REM set calibration_fac = -cf 1.0

@REM set scanner=-s ..\..\config\scanner\HiRezBrainPET_rsmall_Geom_orig
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_orig
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_rsmall_orig.txt

@REM echo ==================================================================================
@REM echo Conversion of HIREZBRAINPET data to CASToR. Scanner: HiRezBrainPET_rsmall_Geom_orig
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope% 

@REM set scanner=-s ..\config\scanner\HiRezBrainPET_rsmall_Geom_z2s
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_z2s
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_rsmall_z2s.txt

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rsmall_Geom_z2s
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope%

@REM set scanner=-s ..\config\scanner\HiRezBrainPET_rsmall_Geom_big
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_big
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_rsmall_big.txt

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rsmall_Geom_big
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope%

@REM set scanner=-s ..\config\scanner\HiRezBrainPET_rsmall_Geom_big_z4s
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_big_z4s
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_rsmall_big_z4s.txt

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rsmall_Geom_big_z4s
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope%

@REM set scanner=-s ..\config\scanner\HiRezBrainPET_rmid_Geom_orig
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_orig
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_rmid_orig.txt

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rmid_Geom_orig
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope%

@REM set scanner=-s ..\config\scanner\HiRezBrainPET_rmid_Geom_z2s
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_z2s
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_rmid_z2s.txt

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rmid_Geom_z2s
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope%

@REM set scanner=-s ..\config\scanner\HiRezBrainPET_rmid_Geom_big_z4s
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_big_z4s
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_rmid_big_z4s.txt

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rmid_Geom_big_z4s
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope%

@REM set scanner=-s ..\config\scanner\HiRezBrainPET_rmid_Geom_big
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_big
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_rmid_big.txt

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rmid_Geom_big
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope%

@REM set scanner=-s ..\config\scanner\HiRezBrainPET_rbig_Geom_orig
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rbig_orig
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_rbig_orig.txt

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rbig_Geom_orig
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope%

@REM set scanner=-s ..\config\scanner\HiRezBrainPET_rbig_Geom_z2s
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rbig_z2s
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_rbig_z2s.txt

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rbig_Geom_z2s
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope%

@REM set scanner=-s ..\config\scanner\HiRezBrainPET_square_Geom
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_square
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_square.txt

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_square_Geom
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope%

@REM set scanner=-s ..\config\scanner\HiRezBrainPET_square_Geom_z2s
@REM set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_square_z2s
@REM set datafile=-il ..\Files\LOR_all_time_ID_Geom_square_z2s.txt

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_square_Geom_z2s
@REM echo ==================================================================================
@REM %command% %datafile% %output% %scanner% %verbose% %isotope%


set scanner=-s ..\config\scanner\HiRezBrainPET_rsmall_Geom_z2s_min90
set output=-o ..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_z2s_min90
set datafile=-il ..\Files\LOR_all_time_ID_Geom_rsmall_z2s.txt

echo ==================================================================================
echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rsmall_Geom_z2s
echo ==================================================================================
%command% %datafile% %output% %scanner% %verbose% %isotope%