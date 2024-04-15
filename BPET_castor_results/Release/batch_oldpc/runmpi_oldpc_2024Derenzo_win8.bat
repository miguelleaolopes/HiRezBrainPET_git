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


set output_1=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_ittest\1test
set output_2=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_ittest\2test
set output_3=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_ittest\3test
set output_4=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_ittest\4test
set output_5=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_ittest\5test
set output_6=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_ittest\6test
set output_7=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_ittest\7test
set output_8=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_ittest\8test
set output_9=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_ittest\9test
set output_10=-dout .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_ittest\10test


set sens=-sens .\..\..\Results\rsmall_z2s\vs015_MLEM_jos_g12_3s\1test\1test_sensitivity.hdr

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
set iteration_1=-it 10:50
@REM set iteration_2=-it 20:25
set iteration_3=-it 25:20
@REM set iteration_4=-it 50:10
@REM set iteration_5=-it 100:5
@REM set iteration_6=-it 250:2
@REM set iteration_7=-it 500:1
set iteration_8=-it 2:50,40:10

set optimizer_MLEM=-opti MLEM
set optimizer_D95=-opti DEPIERRO95
set optimizer_OSL_MLEM=-
set projector=-proj joseph
@REM classicSiddon  distanceDriven  incrementalSiddon  joseph

set penalty_MRF=-pnlt MRF
set penalty_MRP=-pnlt MRP

set penalty_strength_05=-pnlt-beta 0.5
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
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_MLEM_jos_ittest\7test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_7% %iteration_7% %voxels_number_015% %voxels_size_015% %offset% %projector% %psf_3% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM% %sens%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_MLEM_jos_ittest\1test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_1% %iteration_1% %voxels_number_015% %voxels_size_015% %offset% %projector% %psf_3% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM% %sens%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_MLEM_jos_ittest\2test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_2% %iteration_2% %voxels_number_015% %voxels_size_015% %offset% %projector% %psf_3% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM% %sens%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_MLEM_jos_ittest\3test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_3% %iteration_3% %voxels_number_015% %voxels_size_015% %offset% %projector% %psf_3% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM% %sens%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_MLEM_jos_ittest\4test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_4% %iteration_4% %voxels_number_015% %voxels_size_015% %offset% %projector% %psf_3% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM% %sens%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_MLEM_jos_ittest\5test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_5% %iteration_5% %voxels_number_015% %voxels_size_015% %offset% %projector% %psf_3% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM% %sens%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_MLEM_jos_ittest\6test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_6% %iteration_6% %voxels_number_015% 
@REM %voxels_size_015% %offset% %projector% %psf_3% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM% %sens%

@REM echo ==================================================================================
@REM echo Reconstruction with CASToR_Derenzo_all_Geom_rsmall_z2s_df in output rsmall_z2s\vs015_MLEM_jos_ittest\8test
@REM echo ==================================================================================
@REM %mpi_stuff% %recon% %verbose% %datafile2% %output_8% %iteration_8% %voxels_number_015% %voxels_size_015% %offset% %projector% %psf_3% %thread% %last_it% %stats_true% %out_flip% %optimizer_MLEM% %sens%



@echo on
echo ==================================================================================
echo Reconstruction is finished!
echo ==================================================================================