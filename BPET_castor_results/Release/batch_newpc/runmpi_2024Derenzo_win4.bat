@echo off

::  This is the Windows script to run the CASToR project for the HiRezBrainPET Derenzo.
:: To get help about the command-line options, run the program without argument or with '-h', '-help' or '--help' options.

::::::::::::::::::::::::::::
:: Set Command Line Options
::::::::::::::::::::::::::::

set mpi_stuff=mpiexec.exe -n 1
set recon=.\..\castor-recon.exe 

:: And without MPI (comment the one not used)
@REM set recon=.\Release\castor-recon_old.exe

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set verbose=-vb 2
set thread=-th 0
@REM set out_flip=-flip-out Y

set datafile=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_z2s_df.Cdh

set output_1_01=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_g12_3s\1test
set output_1_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_g12_3s\2test

set sens=-sens .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_g12_3s\1test\1test_sensitivity.hdr
set sens2=-sens .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_g12_3s\2test\2test_sensitivity.hdr

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@REM set voxels_number=-dim 300,300,300
set voxels_number=-dim 200,200,40

@REM set fov_size=-fov 300.,300.,300.
@REM set fov_size=-fov 30.,30.,14.

@REM set voxels_size=-fov 1.,1.,1.
set voxels_size=-vox 0.15,0.15,0.35

set offset=-off 0.,-45.,-8.

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@REM set last_it=-oit -1

:: Number of iterations:subsets
set iteration=-it 2:50,40:10
@REM set iteration=-it 1:1,2:50,40:10

set optimizer=-opti MLEM
set optimizer_D95=-opti DEPIERRO95
set optimizer_OSL=-opti OSL

set projector=-proj joseph
set projector_sid=-proj classicSiddon
@REM classicSiddon  distanceDriven  incrementalSiddon  joseph

set penalty_MRF=-pnlt MRF
set penalty_MRP=-pnlt MRP

set penalty_strength_05=-pnlt-beta 0.5
set penalty_strength_01=-pnlt-beta 0.1

set stats_true=-opti-stat

:: image-based PSF using a stationay Gaussian of X1mm transaxial and X2mm axial FWHM with X3 sigmas in
set psf=-conv gaussian,1.2,0.02,3::psf
set psf_2=-conv gaussian,1.2,0.3,3::psf

@REM set post=-conv gaussian,0.4,0.4,3::post

::::::::::::::::::::::::::::
:: Launch the reconstruction
::::::::::::::::::::::::::::

:: Launch the benchmark
echo ==================================================================================
echo Reconstruction is going on. Should take several minutes depending on the hardware.
echo ==================================================================================

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\voxsize02_MLEM_sid_gauss1_3sig\1test
echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile% %output_1_01% %iteration% %voxels_number% %voxels_size% %offset% %projector% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer% %sens%

%mpi_stuff% %recon% %verbose% %datafile% %output_1_02% %iteration% %voxels_number% %voxels_size% %offset% %projector% %psf_2% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer% %sens2%

@echo on
echo ==================================================================================
echo Reconstruction is finished!
echo ==================================================================================