@echo off

::  This is the Windows script to run the CASToR project for the HiRezBrainPET Derenzo.
:: To get help about the command-line options, run the program without argument or with '-h', '-help' or '--help' options.

::::::::::::::::::::::::::::
:: Set Command Line Options
::::::::::::::::::::::::::::

set mpi_stuff=mpiexec.exe -n 1
set recon=.\..\castor-recon.exe 
set recon=.\..\castor-recon_oldpc.exe 

:: And without MPI (comment the one not used)
@REM set recon=.\Release\castor-recon_old.exe

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set verbose=-vb 2
set thread=-th 0
@REM set out_flip=-flip-out Y

set datafile=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_z2s_df.Cdh

set output_1_01=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_sid_g12_3s\1test
set output_1_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\1test
set output_1_03=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_inc_g12_3s\1test

set output_2_01=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_sid_g12_3s\2test
set output_2_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\2test
set output_2_03=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_inc_g12_3s\2test


set output_3_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\3test
set output_4_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\4test
set output_5_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\5test
set output_6_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\6test
set output_7_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\7test
set output_8_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\8test
set output_9_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\9test
set output_10_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\10test
set output_11_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\11test
set output_12_02=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\12test

@REM set sens=-sens .\..\..\Results\rsmall_z2s\vs015_MLEM_sid_g12_3s\1test\1test_sensitivity.hdr
@REM set sens=-sens .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\1test\1test_sensitivity.hdr

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
set projector_dri=-proj distanceDriven
set projector_inc=-proj incrementalSiddon
@REM classicSiddon  distanceDriven  incrementalSiddon  joseph

set penalty_MRF=-pnlt MRF
set penalty_MRP=-pnlt MRP

set penalty_strength_05=-pnlt-beta 0.5
set penalty_strength_01=-pnlt-beta 0.1

set stats_true=-opti-stat

:: image-based PSF using a stationay Gaussian of X1mm transaxial and X2mm axial FWHM with X3 sigmas in
set psf=-conv gaussian,1.2,0.02,3::psf
set psf_2=-conv gaussian,1.2,0.3,3::psf
set psf_3=-conv gaussian,1.,0.3,3::psf
set psf_4=-conv gaussian,0.8,0.3,3::psf
set psf_5=-conv gaussian,1.2,0.3,1.5::psf
set psf_6=-conv gaussian,1.2,0.7,3::psf
set psf_9=-conv gaussian,1.0,0.02,3::psf
set psf_10=-conv gaussian,1.0,0.02,1.5::psf
set psf_11=-conv gaussian,1.2,0.02,1.5::psf
set psf_12=-conv gaussian,1.0,0.3,1.5::psf

set post=-conv gaussian,0.4,0.4,3::post
set post_2=-conv gaussian,1.0,0.35,3::post

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
@REM %mpi_stuff% %recon% %verbose% %datafile% %output_1_01% %iteration% %voxels_number% %voxels_size% %offset% %projector_sid% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_1_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_1_03% %iteration% %voxels_number% %voxels_size% %offset% %projector_inc% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer%


@REM %mpi_stuff% %recon% %verbose% %datafile% %output_2_01% %iteration% %voxels_number% %voxels_size% %offset% %projector_sid% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_2_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_2_03% %iteration% %voxels_number% %voxels_size% %offset% %projector_inc% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer%


@REM %mpi_stuff% %recon% %verbose% %datafile% %output_3_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf% %thread% %post_2% %last_it% %out_flip% %stats_true% %optimizer%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_4_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf_2% %thread% %last_it% %out_flip% %stats_true% %optimizer%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_5_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf_3% %thread% %last_it% %out_flip% %stats_true% %optimizer%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_6_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf_4% %thread% %last_it% %out_flip% %stats_true% %optimizer%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_7_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf_5% %thread% %last_it% %out_flip% %stats_true% %optimizer%

@REM set sens=-sens .\..\..\Results\rsmall_z2s\vs015_MLEM_dri_g12_3s\8test\8test_sensitivity.hdr
@REM %mpi_stuff% %recon% %verbose% %datafile% %output_8_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf_6% %thread% %last_it% %out_flip% %stats_true% %optimizer% %sens%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_9_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf_9% %thread% %last_it% %out_flip% %stats_true% %optimizer%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_10_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf_10% %thread% %last_it% %out_flip% %stats_true% %optimizer%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_11_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf_11% %thread% %last_it% %out_flip% %stats_true% %optimizer%

%mpi_stuff% %recon% %verbose% %datafile% %output_12_02% %iteration% %voxels_number% %voxels_size% %offset% %projector_dri% %psf_12% %thread% %last_it% %out_flip% %stats_true% %optimizer%

@echo on
echo ==================================================================================
echo Reconstruction is finished!
echo ==================================================================================