@echo off

::  This is the Windows script to run the CASToR project for the HiRezBrainPET Derenzo.
:: To get help about the command-line options, run the program without argument or with '-h', '-help' or '--help' options.

::::::::::::::::::::::::::::
:: Set Command Line Options
::::::::::::::::::::::::::::

set mpi_stuff=mpiexec.exe -n 2
set recon=.\..\castor-recon_oldpc.exe 
set recon=.\..\castor-recon.exe 

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

set output_1=-dout .\..\..\Results\rsmall_z2s\vs015_D95_dri_g12_3s\1test
set output_2=-dout .\..\..\Results\rsmall_z2s\vs015_D95_inc_g12_3s\2test
set output_3=-dout .\..\..\Results\rsmall_z2s\vs015_OSL_dri_g12_3s\1test
set output_4=-dout .\..\..\Results\rsmall_z2s\vs015_OSL_dri_g12_3s\2test
set output_5=-dout .\..\..\Results\rsmall_z2s\vs015_OSL_dri_g12_3s\3test
set output_6=-dout .\..\..\Results\rsmall_z2s\vs015_OSL_dri_g12_3s\4test
set output_7=-dout .\..\..\Results\rsmall_z2s\vs015_OSL_dri_g12_3s\5test
set output_8=-dout .\..\..\Results\rsmall_z2s\vs015_OSL_dri_g12_3s\6test
set output_9=-dout .\..\..\Results\rsmall_z2s\vs015_OSL_dri_g12_3s\7test


@REM set sens=-sens .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_g1_3s\1test\1test_sensitivity.hdr
@REM set sens_2=-sens .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_g1_3s\3test\3test_sensitivity.hdr

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
set iteration=-it 2:50,40:10

set optimizer_MLEM=-opti MLEM
set optimizer_D95=-opti DEPIERRO95
set optimizer_OSL=-opti OSL

set projector=-proj joseph
set projector_sid=-proj classicSiddon
set projector_dri=-proj distanceDriven
set projector_inc=-proj incrementalSiddon

@REM classicSiddon  distanceDriven  incrementalSiddon  joseph

set penalty_MRF=-pnlt MRF
set penalty_MRP=-pnlt MRP

set penalty_strength_01=-pnlt-beta 0.1
set penalty_strength_1=-pnlt-beta 1
set penalty_strength_10=-pnlt-beta 10.0
set penalty_strength_100=-pnlt-beta 100.0
set penalty_strength_500=-pnlt-beta 500.0

set stats_true=-opti-stat

:: image-based PSF using a stationay Gaussian of X1mm transaxial and X2mm axial FWHM with X3 sigmas in
set psf=-conv gaussian,1.2,0.02,3::psf

@REM set post=-conv gaussian,0.4,0.4,3::post

::::::::::::::::::::::::::::
:: Launch the reconstruction
::::::::::::::::::::::::::::

:: Launch the benchmark
echo ==================================================================================
echo Reconstruction is going on. Should take several minutes depending on the hardware.
echo ==================================================================================

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\voxsize01_MLEM_jos_g1_3s\2test_rs
echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_1% %iteration% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf% %thread% %post% %last_it% %stats_true% %out_flip% %optimizer_D95% %penalty_MRF% %penalty_strength_100%

@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_2% %iteration% %voxels_number_015% %voxels_size_015% %offset% %projector_inc% %psf% %thread% %post% %last_it% %stats_true% %out_flip% %optimizer_D95% %penalty_MRF% %penalty_strength_100%

@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_3% %iteration% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf% %thread% %post% %last_it% %stats_true% %out_flip% %optimizer_OSL% %penalty_MRF% %penalty_strength_100% 

@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_4% %iteration% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf% %thread% %post% %last_it% %stats_true% %out_flip% %optimizer_OSL% %penalty_MRF% %penalty_strength_10%

@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_7% %iteration% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf% %thread% %post% %last_it% %stats_true% %out_flip% %optimizer_OSL% %penalty_MRF% %penalty_strength_01%

@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_8% %iteration% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf% %thread% %post% %last_it% %stats_true% %out_flip% %optimizer_OSL% %penalty_MRP% %penalty_strength_01%

@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_9% %iteration% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf% %thread% %post% %last_it% %stats_true% %out_flip% %optimizer_OSL% %penalty_MRF% %penalty_strength_1%

@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_5% %iteration% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf% %thread% %post% %last_it% %stats_true% %out_flip% %optimizer_OSL% %penalty_MRP% %penalty_strength_100%

@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_6% %iteration% %voxels_number_015% %voxels_size_015% %offset% %projector_dri% %psf% %thread% %post% %last_it% %stats_true% %out_flip% %optimizer_OSL% %penalty_MRP% %penalty_strength_10%

@echo on
echo ==================================================================================
echo Reconstruction is finished!
echo ==================================================================================