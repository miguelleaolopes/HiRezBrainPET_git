@echo off

::  This is the Windows script to run the CASToR project for the HiRezBrainPET Derenzo.
:: To get help about the command-line options, run the program without argument or with '-h', '-help' or '--help' options.

::::::::::::::::::::::::::::
:: Set Command Line Options
::::::::::::::::::::::::::::

set mpi_stuff=mpiexec.exe -n 20
set recon=.\..\castor-recon.exe 

:: And without MPI (comment the one not used)
@REM set recon=.\Release\castor-recon_old.exe

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set verbose=-vb 2
set thread=-th 20
@REM set out_flip=-flip-out Y

set datafile=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_orig_df.Cdh
set datafile2=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_z2s_df.Cdh
set datafile3=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rsmall_big_df.Cdh
set datafile4=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_orig_df.Cdh
set datafile5=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_z2s_df.Cdh
set datafile6=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_big_z4s_df.Cdh
set datafile7=-df .\..\..\Files\CASToR_data\CASToR_Derenzo_all_Geom_rmid_big_df.Cdh

set output_2_pnl1=-dout .\..\..\Results\rsmall_z2s\vs02_D95_MRF_b05_jos_g1_3s\1test_rs
set output_2_pnl2=-dout .\..\..\Results\rsmall_z2s\vs02_D95_MRF_b05_jos_g1_3s\1test
set output_2_pnl3=-dout .\..\..\Results\rsmall_z2s\vs02_D95_MRF_b01_jos_g1_3s\1test_rs
set output_2_pnl4=-dout .\..\..\Results\rsmall_z2s\vs02_D95_MRF_b01_jos_g1_3s\1test
set output_2_pnl5=-dout .\..\..\Results\rsmall_z2s\vs02_OSL_MRP_b05_jos_g1_3s\1test
set output_2_pnl6=-dout .\..\..\Results\rsmall_z2s\vs02_OSL_MRP_b01_jos_g1_3s\1test

set output_2_3=-dout .\..\..\Results\rsmall_z2s\vs02_MLEM_jos_g1_3s\3test_rs
set output_2_4=-dout .\..\..\Results\rsmall_z2s\vs02_MLEM_jos_g1_3s\4test_rs

set output_2_01=-dout .\..\..\Results\rsmall_z2s\vs01_MLEM_jos_g1_3s\1test


set sens=-sens .\..\..\Results\rsmall_z2s\vs02_MLEM_jos_g1_3s\5test_old\5test_sensitivity.hdr
set sens_2=-sens .\..\..\Results\rsmall_z2s\vs02_D95_MRF_b05_jos_g1_3s\1test\1test_sensitivity.hdr

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@REM set voxels_number=-dim 300,300,300
set voxels_number=-dim 200,200,50
set voxels_number_01=-dim 400,400,50

@REM set fov_size=-fov 300.,300.,300.
@REM set fov_size=-fov 40.,40.,20.

@REM set voxels_size=-fov 1.,1.,1.
set voxels_size=-vox 0.2,0.2,0.4
set voxels_size_01=-vox 0.1,0.1,0.4

set offset=-off 0.,-45.,-8.
set offset_mid=-off 0.,-55.,-8.
set offset_big=-off 0.,-65.,-8.

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@REM set last_it=-oit -1

:: Number of iterations:subsets
@REM set iteration=-it 5:10,5:5
set iteration=-it 2:50,2:10,2:5
set iteration_2=-it 2:10,5:5,2:2,2:1
set iteration_3=-it 10:5,5:2,5:1

set optimizer=-opti MLEM
set optimizer_D95=-opti DEPIERRO95
set optimizer_OSL=-opti OSL

set projector=-proj joseph
@REM classicSiddon  distanceDriven  incrementalSiddon  joseph

set penalty_MRF=-pnlt MRF
set penalty_MRP=-pnlt MRP

set penalty_strength_05=-pnlt-beta 0.5
set penalty_strength_01=-pnlt-beta 0.1

set stats_true=-opti-stat

:: image-based PSF using a stationay Gaussian of X1mm transaxial and X2mm axial FWHM with X3 sigmas in
set psf=-conv gaussian,1,1,3::psf
set psf_2=-conv gaussian,1.5,1.5,3::psf
set psf_3=-conv gaussian,2,2,3::psf

set post=-conv gaussian,0.4,0.4,3::post

::::::::::::::::::::::::::::
:: Launch the reconstruction
::::::::::::::::::::::::::::

:: Launch the benchmark
echo ==================================================================================
echo Reconstruction is going on. Should take several minutes depending on the hardware.
echo ==================================================================================

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_orig_df in output rsmall_orig\vs02_MLEM_jos_gau05_3s\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile% %output% %sens% %iteration% %voxels_number% %voxels_size% %offset% %optimizer% %projector% %psf% %thread% %post% %last_it% %out_flip%

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs02_D95_MRF_b05_jos_g1_3s\1test_rs
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile2% %output_2_pnl1% %iteration% %voxels_number% %voxels_size% %offset% %projector% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer_D95% %penalty_MRF% %penalty_strength_05% %sens%

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs02_D95_MRF_b05_jos_g1_3s\1test
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile2% %output_2_pnl2% %iteration% %voxels_number% %voxels_size% %offset% %projector% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer_D95% %penalty_MRF% %penalty_strength_05%

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs02_D95_MRF_b01_jos_g1_3s\1test_rs
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile2% %output_2_pnl3% %iteration% %voxels_number% %voxels_size% %offset% %projector% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer_D95% %penalty_MRF% %penalty_strength_01% %sens_2%

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs02_D95_MRF_b01_jos_g1_3s\1test
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile2% %output_2_pnl4% %iteration% %voxels_number% %voxels_size% %offset% %projector% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer_D95% %penalty_MRF% %penalty_strength_01%

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs02_OSL_MRP_b05_jos_g1_3s\1test
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile2% %output_2_pnl5% %iteration% %voxels_number% %voxels_size% %offset% %projector% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer_OSL% %penalty_MRP% %penalty_strength_05%

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs02_OSL_MRP_b01_jos_g1_3s\1test
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile2% %output_2_pnl6% %iteration% %voxels_number% %voxels_size% %offset% %projector% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer_OSL% %penalty_MRP% %penalty_strength_01%

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs02_MLEM_jos_g1_3s\3test_rs
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile2% %output_2_3% %iteration_2% %voxels_number% %voxels_size% %offset% %projector% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer% %sens%

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs02_MLEM_jos_g1_3s\4test_rs
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile2% %output_2_4% %iteration_3% %voxels_number% %voxels_size% %offset% %projector% %psf% %thread% %post% %last_it% %out_flip% %stats_true% %optimizer% %sens%

echo ==================================================================================
echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\voxsize01_MLEM_jos_g1_3s\1test
echo ==================================================================================
%mpi_stuff% %recon% %verbose% %datafile2% %output_2_01% %iteration% %voxels_number_01% %voxels_size_01% %offset% %projector% %psf% %thread% %post% %last_it% %stats_true% %out_flip% %optimizer%


@echo on
echo ==================================================================================
echo Reconstruction is finished!
echo ==================================================================================