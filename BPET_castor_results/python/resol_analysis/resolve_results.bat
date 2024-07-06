:: Batch 
@echo off
setlocal enabledelayedexpansion

set compute=python computeDerenzoValleyToPeak.py
set config= -c example/config.json
set method= -m parabola
@REM set shows= --showTriangPos
@REM set shows= --showTriangPos --showSpotsPos --showLinesProfile
set zprofile= -z 20 30

:: File name
set file_arg= -f
set file_path= ..\..\Results\rsmall_z2s\vs015_geoms\
set file_name=1test
:: Add the file path for each file name in the separations
set file_iterations=
for /L %%i in (5,5,50) do (
    set file_iterations=!file_iterations! !file_name!_it%%i.hdr
)
echo %file_iterations%

set file_ends=
for %%i in (%file_iterations%) do (
    echo !file_path!!file_name!\%%i
    set file_ends=!file_ends! !file_path!!file_name!\%%i
)

set output_file=example/output.txt
:: Run the python script

echo Running compute command...
%compute% %config% %method% %shows% %zprofile% %file_arg% %file_ends%
%compute% %config% %method% %shows% %zprofile% %file_arg% %file_ends% > %output_file%

@REM echo on
:: Read the file line by line
for /f "delims=" %%a in (%output_file%) do (
    echo %%a | findstr "Resolvability [%]:" >nul && (
        :: If it does, save the values
        set "average_vpr=%%a"
        echo !average_vpr!
    )
)

endlocal
