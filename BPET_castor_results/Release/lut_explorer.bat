@echo off

set command=.\castor-scannerLUTExplorer.exe
set scanner=-sf ../../BPET_castor_v3.1.1/config/scanner/HiRezBrainPET_rsmall_Geom_orig.hscan
@REM set scanner=-sf ../../BPET_castor_v3.1.1/config/scanner/HiRezBrainPET_rsmall_Geom_z2s.hscan
@REM set scanner=-sf ../../BPET_castor_v3.1.1/config/scanner/HiRezBrainPET_rsmall_Geom_big.hscan
@REM set scanner=-sf ../../BPET_castor_v3.1.1/config/scanner/HiRezBrainPET_rmid_Geom_orig.hscan
@REM set scanner=-sf ../../BPET_castor_v3.1.1/config/scanner/HiRezBrainPET_rmid_Geom_z2s.hscan
@REM set scanner=-sf ../../BPET_castor_v3.1.1/config/scanner/HiRezBrainPET_rmid_Geom_big_z4s.hscan
@REM set scanner=-sf ../../BPET_castor_v3.1.1/config/scanner/HiRezBrainPET_rmid_Geom_big.hscan
@REM set scanner=-sf ../../BPET_castor_v3.1.1/config/scanner/HiRezBrainPET_rbig_Geom_orig.hscan
@REM set scanner=-sf ../../BPET_castor_v3.1.1/config/scanner/HiRezBrainPET_rbig_Geom_z2s.hscan
@REM set scanner=-sf ../../BPET_castor_v3.1.1/config/scanner/HiRezBrainPET_square_Geom.hscan
@REM set scanner=-sf ../../BPET_castor_v3.1.1/config/scanner/HiRezBrainPET_sqaure_Geom_z2s.hscan

%command% %scanner% -e
