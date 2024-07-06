@echo off

::  This is the Windows script to run the CASToR project for the HiRezBrainPET Derenzo.
:: To get help about the command-line options, run the program without argument or with '-h', '-help' or '--help' options.

::::::::::::::::::::::::::::
:: Set Command Line Options
::::::::::::::::::::::::::::

set mpi_stuff=mpiexec.exe
set recon=.\..\castor-recon_oldpc.exe 
set recon=.\..\castor-recon.exe 

:: And without MPI (comment the one not used)
@REM set recon=.\Release\castor-recon_old.exe

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set verbose=-vb 2
set thread=-th 2
@REM set out_flip=-flip-out Y

set datafile2=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_z2s_df.Cdh
@REM set datafile2=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_z2s_min90_df.Cdh

set output_1=-dout .\..\..\Results\rsmall_z2s\vs015_min90\opt

set sens=-sens .\..\..\Results\rsmall_z2s\vs015_min90\opt\opt_sensitivity.hdr

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
set iteration_4=-it 50:10
@REM set iteration_8=-it 2:50,40:10

set optimizer_MLEM=-opti MLEM
set projector=-proj joseph
set projector_dri=-proj distanceDriven
set projector_inc=-proj incrementalSiddon

set penalty_MRF=-pnlt MRF
set penalty_MRP=-pnlt MRP

set penalty_strength_03=-pnlt-beta 0.3
set penalty_strength_05=-pnlt-beta 0.5
set penalty_strength_08=-pnlt-beta 0.8
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

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_min90\opt
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile2% %output_1% %iteration_4% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf_3% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM%
