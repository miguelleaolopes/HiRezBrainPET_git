@echo off

@REM This script converts the txt data with the positions of the LOR (t,x1,y1,z1,x2,y2,z2) to crystal ID.
@REM The data is stored in the file LOR_all_points_time.txt



@REM Define variables for the command line
set command=.\castor-txtConversionCrystalsID_BPET.exe
set verbose=-vb 3
set txtfile=-txt ..\Files\LOR_all_points_time.txt

set output=-o LOR_all_time_ID_Geom_rsmall_orig.txt 
set scanner=-sf ..\config\scanner\HiRezBrainPET_rsmall_Geom_orig

echo ==================================================================================
echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rsmall_Geom_orig
echo ==================================================================================
@REM %command% %txtfile% %output% %scanner% %verbose%

@REM set output=-o LOR_all_time_ID_Geom_rsmall_z2s.txt
@REM set scanner=-sf ..\config\scanner\HiRezBrainPET_rsmall_Geom_z2s

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rsmall_Geom_z2s
@REM echo ==================================================================================
@REM %command% %txtfile% %output% %scanner% %verbose%

@REM set output=-o LOR_all_time_ID_Geom_rsmall_big.txt
@REM set scanner=-sf ..\config\scanner\HiRezBrainPET_rsmall_Geom_big

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rsmall_Geom_big
@REM echo ==================================================================================
@REM %command% %txtfile% %output% %scanner% %verbose%

@REM set output=-o LOR_all_time_ID_Geom_rmid_orig.txt
@REM set scanner=-sf ..\config\scanner\HiRezBrainPET_rmid_Geom_orig

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rmid_Geom_orig
@REM echo ==================================================================================
@REM %command% %txtfile% %output% %scanner% %verbose%

@REM set output=-o LOR_all_time_ID_Geom_rmid_z2s.txt
@REM set scanner=-sf ..\config\scanner\HiRezBrainPET_rmid_Geom_z2s

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rmid_Geom_z2s
@REM echo ==================================================================================
@REM %command% %txtfile% %output% %scanner% %verbose%

@REM set output=-o LOR_all_time_ID_Geom_rmid_big_z4s.txt
@REM set scanner=-sf ..\config\scanner\HiRezBrainPET_rmid_Geom_big_z4s

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rmid_Geom_big_z4s
@REM echo ==================================================================================
@REM %command% %txtfile% %output% %scanner% %verbose%

@REM set output=-o LOR_all_time_ID_Geom_rmid_big.txt
@REM set scanner=-sf ..\config\scanner\HiRezBrainPET_rmid_Geom_big

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rmid_Geom_big
@REM echo ==================================================================================
@REM %command% %txtfile% %output% %scanner% %verbose%

@REM set output=-o LOR_all_time_ID_Geom_rbig_orig.txt
@REM set scanner=-sf ..\config\scanner\HiRezBrainPET_rbig_Geom_orig

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rbig_Geom_orig
@REM echo ==================================================================================
@REM %command% %txtfile% %output% %scanner% %verbose%

@REM set output=-o LOR_all_time_ID_Geom_rbig_z2s.txt
@REM set scanner=-sf ..\config\scanner\HiRezBrainPET_rbig_Geom_z2s

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_rbig_Geom_z2s
@REM echo ==================================================================================
@REM %command% %txtfile% %output% %scanner% %verbose%

@REM set output=-o LOR_all_time_ID_Geom_square.txt
@REM set scanner=-sf ..\config\scanner\HiRezBrainPET_square_Geom

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_square_Geom
@REM echo ==================================================================================
@REM %command% %txtfile% %output% %scanner% %verbose%

@REM set output=-o LOR_all_time_ID_Geom_square_z2s.txt
@REM set scanner=-sf ..\config\scanner\HiRezBrainPET_square_Geom_z2s

@REM echo ==================================================================================
@REM echo Conversion of TXT file with crystal ID. Scanner: HiRezBrainPET_square_Geom_z2s
@REM echo ==================================================================================
@REM %command% %txtfile% %output% %scanner% %verbose%
