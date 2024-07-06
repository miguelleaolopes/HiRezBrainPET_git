@echo off

::  This is the Windows script to run the CASToR project for the HiRezBrainPET Derenzo.
:: To get help about the command-line options, run the program without argument or with '-h', '-help' or '--help' options.

::::::::::::::::::::::::::::
:: Set Command Line Options
::::::::::::::::::::::::::::

set mpi_stuff=mpiexec.exe
set recon=.\..\castor-recon.exe 
set recon=.\..\castor-recon_oldpc.exe 

:: And without MPI (comment the one not used)
@REM set recon=.\Release\castor-recon_old.exe

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set verbose=-vb 2
set thread=-th 0
@REM set out_flip=-flip-out Y

set datafile=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_orig_df.Cdh
set datafile2=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_z2s_df.Cdh
set datafile3=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_big_df.Cdh
set datafile4=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_orig_df.Cdh
set datafile5=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_z2s_df.Cdh
set datafile6=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_big_z4s_df.Cdh
set datafile7=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_big_df.Cdh


set output_1=-dout .\..\..\Results\rsmall_z2s\vs015_opts_dri\MLEM
set output_2=-dout .\..\..\Results\rsmall_z2s\vs015_opts_dri\D95\b1
set output_2_2=-dout .\..\..\Results\rsmall_z2s\vs015_opts_dri\D95\b1_rs
set output_3=-dout .\..\..\Results\rsmall_z2s\vs015_opts_dri\OSL_P\b1
set output_4=-dout .\..\..\Results\rsmall_z2s\vs015_opts_dri\OSL_F\b1



set sens=-sens .\..\..\Results\rsmall_z2s\vs015_MLEM_projs\dri\dri_sensitivity.hdr

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set voxels_number_015=-dim 200,200,40

@REM set fov_size_02=-fov 30.,30.,14.

set voxels_size_015=-vox 0.15,0.15,0.35

set offset=-off 0.,-45.,-8.
set offset_mid=-off 0.,-55.,-8.
set offset_big=-off 0.,-65.,-8.

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@REM set last_it=-oit -1

:: Number of iterations:subsets
@REM set iteration=-it 5:10,5:5
@REM set iteration_1=-it 10:50
@REM set iteration_2=-it 20:25
@REM set iteration_3=-it 25:20
@REM set iteration_4=-it 50:10
@REM set iteration_5=-it 100:5
@REM set iteration_6=-it 250:2
@REM set iteration_7=-it 500:1
set iteration_8=-it 2:50,40:10

set optimizer_MLEM=-opti MLEM
set optimizer_D95=-opti DEPIERRO95
set optimizer_OSL=-opti OSL
set projector=-proj joseph
set projector_sid=-proj classicSiddon
set projector_dri=-proj distanceDriven
set projector_inc=-proj incrementalSiddon

set penalty_MRF=-pnlt MRF
set penalty_MRP=-pnlt MRP

set penalty_strength_05=-pnlt-beta 0.5
set penalty_strength_1=-pnlt-beta 1
set penalty_strength_10=-pnlt-beta 10.0
set penalty_strength_100=-pnlt-beta 100.0
set penalty_strength_500=-pnlt-beta 500.0
set penalty_strength_01=-pnlt-beta 0.1

set stats_true=-opti-stat

:: image-based PSF using a stationay Gaussian of X1mm transaxial and X2mm axial FWHM with X3 sigmas in
set psf=-conv gaussian,1,1,3::psf
set psf_2=-conv gaussian,1,1,2::psf
set psf_3=-conv gaussian,1.2,0.02,3::psf

set post=-conv gaussian,0.4,0.4,3::post

::::::::::::::::::::::::::::
:: Launch the reconstruction
::::::::::::::::::::::::::::

:: Launch the benchmark
echo ==================================================================================
echo Reconstruction is going on. Should take several minutes depending on the hardware.
echo ==================================================================================

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_opts_dri\MLEM
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_1% %iteration_8% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf_3% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM% %sens%

@REM set sensb1=-sens .\..\..\Results\rsmall_z2s\vs015_opts_dri\D95\b1\b1_sensitivity.hdr
@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_opts_dri\D95\b1
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_2% %iteration_8% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf_3% %thread% %stats_true% %out_flip% %optimizer_D95% %penalty_MRF% %penalty_strength_1% %sensb1%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_opts_dri\OSL_P\b1
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_3% %iteration_8% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf_3% %thread% %stats_true% %out_flip% %optimizer_OSL% %penalty_MRP% %penalty_strength_1%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_opts_dri\OSL_F\b1
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_4% %iteration_8% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf_3% %thread% %stats_true% %out_flip% %optimizer_OSL% %penalty_MRF% %penalty_strength_1%

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_opts_dri\D95\b1_rs
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile2% %output_2_2% %iteration_8% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf_3% %thread% %stats_true% %out_flip% %optimizer_D95% %penalty_MRF% %penalty_strength_1% %sens%

@echo on
echo ==================================================================================
echo Reconstruction is finished!
echo ==================================================================================