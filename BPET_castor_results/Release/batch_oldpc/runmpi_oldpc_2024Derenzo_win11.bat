@echo off

::  This is the Windows script to run the CASToR project for the HiRezBrainPET Derenzo.
:: To get help about the command-line options, run the program without argument or with '-h', '-help' or '--help' options.

::::::::::::::::::::::::::::
:: Set Command Line Options
::::::::::::::::::::::::::::

set mpi_stuff=mpiexec.exe -n 2
set recon=.\..\castor-recon.exe 
set recon=.\..\castor-recon_oldpc.exe 

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set verbose=-vb 2
set thread=-th 0
@REM set out_flip=-flip-out Y

set datafile1=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_orig_df.Cdh
set datafile2=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_z2s_df.Cdh
set datafile3=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_big_df.Cdh
set datafile4=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_orig_df.Cdh
set datafile5=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_z2s_df.Cdh
set datafile6=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_big_z4s_df.Cdh
set datafile7=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_big_df.Cdh
set datafile8=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rbig_orig_df.Cdh
set datafile9=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rbig_z2s_df.Cdh
set datafile10=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_big_z4s_df.Cdh


set output1=-dout .\..\..\Results\rsmall_orig\vs015_geoms\1test
set output2=-dout .\..\..\Results\rsmall_z2s\vs015_geoms\1test
set output3=-dout .\..\..\Results\rsmall_big\vs015_geoms\1test
set output4=-dout .\..\..\Results\rmid_orig\vs015_geoms\1test
set output5=-dout .\..\..\Results\rmid_z2s\vs015_geoms\1test
set output6=-dout .\..\..\Results\rmid_big_z4s\vs015_geoms\1test
set output7=-dout .\..\..\Results\rmid_big\vs015_geoms\1test
set output8=-dout .\..\..\Results\rbig_orig\vs015_geoms\1test
set output9=-dout .\..\..\Results\rbig_z2s\vs015_geoms\1test
set output10=-dout .\..\..\Results\rsmall_big_z4s\vs015_geoms\1test


set sens2=-sens .\..\..\Results\rsmall_z2s\vs015_geoms\1test\1test_sensitivity.hdr

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set voxels_number_015=-dim 200,200,40

@REM set fov_size_02=-fov 30.,30.,14.

set voxels_size_015=-vox 0.15,0.15,0.35

@REM TEMPORARY
set voxels_number_015=-dim 220,220,45
set voxels_size_015=-vox 0.15,0.15,0.35

set offset_small=-off 0.,-45.,-8.
set offset_mid=-off 0.,-45.,-8.
set offset_big=-off 0.,-45.,-8.

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@REM set last_it=-oit -1

:: Number of iterations:subsets
set iteration=-it 50:10
@REM set iteration=-it 2:50,40:10

set optimizer_MLEM=-opti MLEM

set projector_dri=-proj distanceDriven

set stats_true=-opti-stat

:: image-based PSF using a stationay Gaussian of X1mm transaxial and X2mm axial FWHM with X3 sigmas in
set psf=-conv gaussian,1.0,0.02,3::psf

@REM set post=-conv gaussian,0.4,0.4,3::post

::::::::::::::::::::::::::::
:: Launch the reconstruction
::::::::::::::::::::::::::::

:: Launch the benchmark
echo ==================================================================================
echo Reconstruction is going on. Should take several minutes depending on the hardware.
echo ==================================================================================

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_orig_df in output rsmall_orig\vs015_geoms\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile1% %output1% %iteration% %voxels_number_015% %voxels_size_015% %offset_small% %projector_dri% %psf% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_geoms\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output2% %iteration% %voxels_number_015% %voxels_size_015% %offset_small% %projector_dri% %psf% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_big_df in output rsmall_big\vs015_geoms\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile3% %output3% %iteration% %voxels_number_015% %voxels_size_015% %offset_small% %projector_dri% %psf% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_big_z4s_df in output rsmall_big_z4s\vs015_geoms\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile10% %output10% %iteration% %voxels_number_015% %voxels_size_015% %offset_small% %projector_dri% %psf% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_orig_df in output rmid_orig\vs015_geoms\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile4% %output4% %iteration% %voxels_number_015% %voxels_size_015% %offset_mid% %projector_dri% %psf% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_z2s_df in output rmid_z2s\vs015_geoms\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile5% %output5% %iteration% %voxels_number_015% %voxels_size_015% %offset_mid% %projector_dri% %psf% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM%

@REM set sens=-sens .\..\..\Results\rmid_big_z4s\vs015_geoms\1test\1test_sensitivity.hdr
@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_big_z4s_df in output rmid_big_z4s\vs015_geoms\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile6% %output6% %iteration% %voxels_number_015% %voxels_size_015% %offset_mid% %projector_dri% %psf% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM% %sens%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rbig_orig_df in output rbig_orig\vs015_geoms\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile8% %output8% %iteration% %voxels_number_015% %voxels_size_015% %offset_big% %projector_dri% %psf% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rbig_z2s_df in output rbig_z2s\vs015_geoms\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile9% %output9% %iteration% %voxels_number_015% %voxels_size_015% %offset_big% %projector_dri% %psf% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM%

set sens=-sens .\..\..\Results\rmid_big\vs015_geoms\1test\1test_sensitivity.hdr
echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_big_df in output rmid_big\vs015_geoms\1test
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile7% %output7% %iteration% %voxels_number_015% %voxels_size_015% %offset_mid% %projector_dri% %psf% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM% %sens%
