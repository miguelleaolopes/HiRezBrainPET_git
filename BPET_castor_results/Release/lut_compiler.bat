@echo off

@REM set command= .\castor-PetScannerLutEx.exe	
set command= .\..\..\BPET_castor_v3.1.1\build\Release\castor-PetScannerLutEx_BPET.exe	

:: Constant parameters
set defaultfov= -defaultfov 300,300
set minangle= -minangle 22.5
set startangle= -startangle 90

:: HiRezBrainPET_rsmall_Geom_orig
set alias= -alias HiRezBrainPET_rsmall_Geom_orig
set radius= -radius 150 
set nbrsectors= -nbrsectors 942
set nbcrystals= -nbcrystals 1,300
set sizecrystals= -sizecrystals 10,1,1
set defaultdim= -defaultdim 300,300

@REM :: HiRezBrainPET_rsmall_Geom_z2s
@REM set alias= -alias HiRezBrainPET_rsmall_Geom_z2s
@REM set radius= -radius 150 
@REM set nbrsectors= -nbrsectors 1884
@REM set nbcrystals= -nbcrystals 1,150
@REM set sizecrystals= -sizecrystals 10,0.5,2
@REM set defaultdim= -defaultdim 600,150

@REM :: HiRezBrainPET_rsmall_Geom_big
@REM set alias= -alias HiRezBrainPET_rsmall_Geom_big
@REM set radius= -radius 150 
@REM set nbrsectors= -nbrsectors 3768
@REM set nbcrystals= -nbcrystals 1,150
@REM set sizecrystals= -sizecrystals 10,0.25,2
@REM set defaultdim= -defaultdim 1200,150

@REM :: HiRezBrainPET_rmid_Geom_orig
@REM set alias= -alias HiRezBrainPET_rmid_Geom_orig
@REM set radius= -radius 181.07
@REM set nbrsectors= -nbrsectors 1138
@REM set nbcrystals= -nbcrystals 1,300
@REM set sizecrystals= -sizecrystals 10,1,1
@REM set defaultdim= -defaultdim 300,300

@REM :: HiRezBrainPET_rmid_Geom_z2s
@REM set alias= -alias HiRezBrainPET_rmid_Geom_z2s
@REM set radius= -radius 181.07
@REM set nbrsectors= -nbrsectors 2276
@REM set nbcrystals= -nbcrystals 1,150
@REM set sizecrystals= -sizecrystals 10,0.5,2
@REM set defaultdim= -defaultdim 600,150

@REM :: HiRezBrainPET_rmid_Geom_big_z4s
@REM set alias= -alias HiRezBrainPET_rmid_Geom_big_z4s
@REM set radius= -radius 181.07
@REM set nbrsectors= -nbrsectors 4552
@REM set nbcrystals= -nbcrystals 1,75
@REM set sizecrystals= -sizecrystals 10,0.25,4
@REM set defaultdim= -defaultdim 1200,75

@REM :: HiRezBrainPET_rmid_Geom_big
@REM set alias= -alias HiRezBrainPET_rmid_Geom_big
@REM set radius= -radius 181.07
@REM set nbrsectors= -nbrsectors 4552
@REM set nbcrystals= -nbcrystals 1,150
@REM set sizecrystals= -sizecrystals 10,0.25,2
@REM set defaultdim= -defaultdim 1200,150

@REM :: HiRezBrainPET_rbig_Geom_orig
@REM set alias= -alias HiRezBrainPET_rbig_Geom_orig
@REM set radius= -radius 212.13
@REM set nbrsectors= -nbrsectors 1332
@REM set nbcrystals= -nbcrystals 1,300
@REM set sizecrystals= -sizecrystals 10,1,1
@REM set defaultdim= -defaultdim 300,300

@REM :: HiRezBrainPET_rbig_Geom_z2s
@REM set alias= -alias HiRezBrainPET_rbig_Geom_z2s
@REM set radius= -radius 212.13
@REM set nbrsectors= -nbrsectors 2664
@REM set nbcrystals= -nbcrystals 1,150
@REM set sizecrystals= -sizecrystals 10,0.5,2
@REM set defaultdim= -defaultdim 600,150


@REM :: HiRezBrainPET_square_Geom
@REM set alias= -alias HiRezBrainPET_square_Geom
@REM set radius= -radius 145
@REM set nbrsectors= -nbrsectors 4
@REM set nbcrystals= -nbcrystals 300,300
@REM set sizecrystals= -sizecrystals 10,1,1
@REM set defaultdim= -defaultdim 300,300
@REM set startangle= -startangle 45

:: HiRezBrainPET_square_Geom_z2s
@REM set alias= -alias HiRezBrainPET_square_Geom_z2s
@REM set radius= -radius 145
@REM set nbrsectors= -nbrsectors 4
@REM set nbcrystals= -nbcrystals 600,150
@REM set sizecrystals= -sizecrystals 10,0.5,2
@REM set defaultdim= -defaultdim 600,150
@REM set startangle= -startangle 45

echo ==================================================================================
echo LUT Compiler
echo ==================================================================================

%command% %alias% %radius% %nbrsectors% %nbcrystals% %sizecrystals% %defaultdim% %defaultfov% %minangle% %startangle%
@REM can you change so that it computes the scanner elements position for the first rsector, not specific to the right side of isocenter, but in general given the rsector_first_angle
