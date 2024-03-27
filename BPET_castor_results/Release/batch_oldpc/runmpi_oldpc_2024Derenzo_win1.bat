@echo off

::  This is the Windows script to run the CASToR project for the HiRezBrainPET Derenzo.
:: To get help about the command-line options, run the program without argument or with '-h', '-help' or '--help' options.

::::::::::::::::::::::::::::
:: Set Command Line Options
::::::::::::::::::::::::::::

set mpi_stuff=mpiexec.exe
@REM -n 4
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
set datafile8=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rbig_orig_df.Cdh
set datafile9=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rbig_z2s_df.Cdh

set output=-dout .\..\..\Results\rsmall_orig\vs02_MLEM_jos_gauss05_3sig\1test
set output2=-dout .\..\..\Results\rsmall_z2s\vs02_MLEM_jos_gauss05_3sig\1test
set output3=-dout .\..\..\Results\rsmall_big\vs02_MLEM_jos_gauss05_3sig\1test
set output4=-dout .\..\..\Results\rmid_orig\vs02_MLEM_jos_gauss05_3sig\1test
set output5=-dout .\..\..\Results\rmid_z2s\vs02_MLEM_jos_gauss05_3sig\1test
set output6=-dout .\..\..\Results\rmid_z4s\vs02_MLEM_jos_gauss05_3sig\1test
set output7=-dout .\..\..\Results\rmid_big\vs02_MLEM_jos_gauss05_3sig\1test
set output8=-dout .\..\..\Results\rbig_orig\vs02_MLEM_jos_gauss05_3sig\1test
set output9=-dout .\..\..\Results\rbig_z2s\vs02_MLEM_jos_gauss05_3sig\1test

@REM set output_2=-dout .\..\..\Results\1test\all_z2s_fov02_OSL_1t
@REM set output_3=-dout .\..\..\Results\1test\all_z2s_fov02_DEPIERRO95_1t

@REM set sens=-sens .\..\..\Results\1test\all_z2s_fov02_MLEM_1t\all_z2s_fov02_MLEM_1t_sensitivity.hdr

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@REM set voxels_number=-dim 300,300,300
set voxels_number=-dim 200,200,50

@REM set fov_size=-fov 300.,300.,300.
@REM set fov_size=-fov 40.,40.,20.

@REM set voxels_size=-fov 1.,1.,1.
set voxels_size=-vox 0.2,0.2,0.4

set offset=-off 0.,-45.,-8.
set offset_mid=-off 0.,-55.,-8.
set offset_big=-off 0.,-65.,-8.`

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@REM set last_it=-oit -1

:: Number of iterations:subsets
set iteration=-it 5:10,5:5
set iteration=-it 2:50,2:10

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

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_orig_df in output rsmall_orig\vs02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile% %output% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output2% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_big_df in output rsmall_big\vs02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile3% %output3% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_orig_df in output rmid_orig\vs02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile4% %output4% %sens% %iteration% %voxels_number% %voxels_size% %offset_mid% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_z2s_df in output rmid_z2s\vs02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile5% %output5% %sens% %iteration% %voxels_number% %voxels_size% %offset_mid% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_big_z4s_df in output rmid_z4s\vs02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile6% %output6% %sens% %iteration% %voxels_number% %voxels_size% %offset_mid% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rmid_big_df in output rmid_big\vs02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile7% %output7% %sens% %iteration% %voxels_number% %voxels_size% %offset_mid% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rbig_orig_df in output rbig_orig\vs02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile8% %output8% %sens% %iteration% %voxels_number% %voxels_size% %offset_big% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rbig_z2s_df in output rbig_z2s\vs02_MLEM_jos_gauss05_3sig\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile9% %output9% %sens% %iteration% %voxels_number% %voxels_size% %offset_big% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%


@echo on
echo ==================================================================================
echo Reconstruction is finished!
echo ==================================================================================