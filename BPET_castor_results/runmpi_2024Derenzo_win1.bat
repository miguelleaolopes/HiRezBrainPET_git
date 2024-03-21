@echo off

::  This is the Windows script to run the CASToR project for the HiRezBrainPET Derenzo.
:: To get help about the command-line options, run the program without argument or with '-h', '-help' or '--help' options.

::::::::::::::::::::::::::::
:: Set Command Line Options
::::::::::::::::::::::::::::

set mpi_stuff=mpiexec.exe
@REM -n 15
set recon=.\Release\castor-recon.exe 

:: And without MPI (comment the one not used)
@REM set recon=.\Release\castor-recon_old.exe

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set verbose=-vb 2
set thread=-th 0
@REM set out_flip=-flip-out Y

set datafile=-df Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_orig_df.Cdh
@REM set datafile2=-df Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_z2s_df.Cdh
@REM set datafile3=-df Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_big_df.Cdh
@REM set datafile4=-df Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_orig_df.Cdh
@REM set datafile5=-df Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_z2s_df.Cdh
@REM set datafile6=-df Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_big_z4s_df.Cdh
@REM set datafile7=-df Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_big_df.Cdh

set output=-dout Results\2024_Derenzo\rsmall_orig\voxsize02_MLEM_jos_gauss05_3sig_1t
@REM set output2=-dout Results\2024_Derenzo\rsmall_z2s\voxsize02_MLEM_jos_gauss05_3sig\1test
@REM set output3=-dout Results\2024_Derenzo\rsmall_big\voxsize02_MLEM_jos_gauss05_3sig\1test
@REM set output4=-dout Results\2024_Derenzo\rmid_orig\voxsize02_MLEM_jos_gauss05_3sig\1test
@REM set output5=-dout Results\2024_Derenzo\rmid_z2s\voxsize02_MLEM_jos_gauss05_3sig\1test
@REM set output6=-dout Results\2024_Derenzo\rmid_z4s\voxsize02_MLEM_jos_gauss05_3sig\1test
@REM set output7=-dout Results\2024_Derenzo\rmid_big\voxsize02_MLEM_jos_gauss05_3sig\1test

@REM set output_2=-douResults\2024_Derenzo\1test\all_z2s_fov02_OSL_1t
@REM set output_3=-douResults\2024_Derenzo\1test\all_z2s_fov02_DEPIERRO95_1t

@REM set sens=-senResults\2024_Derenzo\1test\all_z2s_fov02_MLEM_1t\all_z2s_fov02_MLEM_1t_sensitivity.hdr

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@REM set last_it=-oit -1

:: Number of iterations:subsets
set iteration=-it 5:10,5:5
set iteration=-it 2:50,2:10

@REM set voxels_number=-dim 300,300,300
set voxels_number=-dim 200,200,50

@REM set fov_size=-fov 300.,300.,300.
@REM set fov_size=-fov 40.,40.,20.

@REM set voxels_size=-fov 1.,1.,1.
set voxels_size=-vox 0.2,0.2,0.4

set offset=-off 0.,-45.,-8.

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


set optimizer=-opti MLEM
set optimizer_2=-opti OSL
set optimizer_3=-opti DEPIERRO95
@REM OSL  MLEM   DEPIERRO95

set projector=-proj joseph
@REM classicSiddon  distanceDriven  incrementalSiddon  joseph


:: image-based PSF using a stationay Gaussian of X1mm transaxial and X2mm axial FWHM with X3 sigmas in
set psf=-conv gaussian,0.5,0.5,3::psf
set psf_2=-conv gaussian,0.5,0.5,3::psf
set psf_3=-conv gaussian,2.,2.,3::psf

set post=-conv gaussian,0.4,0.4,3::post

::::::::::::::::::::::::::::
:: Launch the reconstruction
::::::::::::::::::::::::::::

:: Launch the benchmark
echo ==================================================================================
echo Reconstruction is going on. Should take several minutes depending on the hardware.
echo ==================================================================================

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_orig_df in output rsmall_orig\voxsize02_MLEM_jos_gauss05_3sig\1test
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile% %output% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\voxsize02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output2% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_big_df in output rsmall_big\voxsize02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile3% %output3% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_orig_df in output rmid_orig\voxsize02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile4% %output4% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_z2s_df in output rmid_z2s\voxsize02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile5% %output5% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_big_z4s_df in output rmid_z4s\voxsize02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile6% %output6% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_big_df in output rmid_big\voxsize02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile7% %output7% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%



@REM %mpi_stuff% %recon% %verbose% %datafile% %output_2% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer_2% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM %mpi_stuff% %recon% %verbose% %datafile% %output_3% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer_3% %projector% %psf% %thread% %post% %last_it% %out_flip%

@echo on
echo ==================================================================================
echo Reconstruction is finished!
echo ==================================================================================