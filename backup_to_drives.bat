@echo off
setlocal enabledelayedexpansion

rem Merges D:\Charlotte source folders into the F: and G: backup drives.
rem
rem Behavior (per folder, per destination drive):
rem   - Files that don't exist on the destination are copied.
rem   - Files that exist on both are only overwritten if the D: copy is newer.
rem   - Nothing is ever deleted from the destination (safe incremental merge).
rem   - Destination folders that don't already exist on F:/G: are created and
rem     flagged as NEW, so you can spot name mismatches (e.g. D:'s "N01-015"
rem     vs an existing "N01-015 Template" on F:/G:) before trusting the merge.
rem
rem Usage:
rem   backup_to_drives.bat            (real run)
rem   backup_to_drives.bat -whatif    (dry run - lists what would copy, copies nothing)

set "SOURCE1=D:\Charlotte\ABAQUS"
set "SOURCE2=D:\Charlotte\cmb247 - Segmentation"
set "LOOSEROOT=D:\"
set "DEST1=F:\"
set "DEST2=G:\"

rem Junk/system files to skip when copying loose files sitting at D:\ root.
set "LOOSE_EXCLUDE=*.dmp *.dmp.xml DumpStack.log.tmp ~$*.xlsx backup_to_drives.bat *.zip"

set "WHATIF="
if /I "%~1"=="-whatif" set "WHATIF=/L"

set "LOGDIR=D:\backup_logs"
if not exist "%LOGDIR%" mkdir "%LOGDIR%"
set "LOGFILE=%LOGDIR%\backup_log.txt"

echo ============================================== >> "%LOGFILE%"
echo Backup run started: %date% %time% >> "%LOGFILE%"
if defined WHATIF echo (DRY RUN - no files will be copied) >> "%LOGFILE%"
echo ============================================== >> "%LOGFILE%"

set "ERRORCOUNT=0"
set "NEWCOUNT=0"

for %%S in ("%SOURCE1%" "%SOURCE2%") do (
    call :ProcessSource "%%~S"
)

call :ProcessLooseFiles

goto :Summary

:ProcessSource
set "SRCROOT=%~1"
if not exist "%SRCROOT%" (
    echo ERROR: source folder not found: %SRCROOT%
    echo ERROR: source folder not found: %SRCROOT% >> "%LOGFILE%"
    goto :eof
)
for /d %%F in ("%SRCROOT%\*") do (
    call :ProcessFolder "%%~fF" "%%~nxF"
)
goto :eof

:ProcessFolder
set "SRC=%~1"
set "NAME=%~2"
set "MISSINGON="

for %%D in ("%DEST1%" "%DEST2%") do (
    if not exist "%%~D%NAME%" set "MISSINGON=!MISSINGON! %%~D"
)

if defined MISSINGON (
    set /a NEWCOUNT+=1
    echo WARNING: NEW FOLDER - "%NAME%" does not exist yet on:!MISSINGON! Check for name mismatches ^(e.g. N01-015 vs "N01-015 Template"^) before trusting this as a merge.
    echo WARNING: NEW FOLDER - "%NAME%" does not exist yet on:!MISSINGON! Check for name mismatches. >> "%LOGFILE%"
)

for %%D in ("%DEST1%" "%DEST2%") do (
    set "DEST=%%~D%NAME%"

    echo Copying "%SRC%" -^> "!DEST!" ...
    echo Copying "%SRC%" -^> "!DEST!" ... >> "%LOGFILE%"
    robocopy "%SRC%" "!DEST!" /E /XO /R:2 /W:5 /NP /NDL /NFL /LOG+:"%LOGFILE%" !WHATIF!
    set "RC=!ERRORLEVEL!"

    if !RC! GEQ 8 (
        set /a ERRORCOUNT+=1
        echo ERROR: robocopy failed for "!DEST!" with exit code !RC!
        echo ERROR: robocopy failed for "!DEST!" with exit code !RC! >> "%LOGFILE%"
    )
)
goto :eof

:ProcessLooseFiles
echo Copying loose files at "%LOOSEROOT%" ...
echo Copying loose files at "%LOOSEROOT%" ... >> "%LOGFILE%"

for %%D in ("%DEST1%" "%DEST2%") do (
    echo   -^> %%~D
    echo   -^> %%~D >> "%LOGFILE%"
    robocopy %LOOSEROOT% %%~D /XF %LOOSE_EXCLUDE% /XO /R:2 /W:5 /NP /NDL /NFL /LOG+:"%LOGFILE%" !WHATIF!
    set "RC=!ERRORLEVEL!"

    if !RC! GEQ 8 (
        set /a ERRORCOUNT+=1
        echo ERROR: robocopy failed for loose files -^> "%%~D" with exit code !RC!
        echo ERROR: robocopy failed for loose files -^> "%%~D" with exit code !RC! >> "%LOGFILE%"
    )
)
goto :eof

:Summary
echo.
echo ----- Summary -----
echo New folders created on a destination: %NEWCOUNT%
echo Folders with errors: %ERRORCOUNT%
echo Full log: %LOGFILE%
echo ----- Summary ----- >> "%LOGFILE%"
echo New folders created on a destination: %NEWCOUNT% >> "%LOGFILE%"
echo Folders with errors: %ERRORCOUNT% >> "%LOGFILE%"
echo Backup run finished: %date% %time% >> "%LOGFILE%"

:End
echo.
pause
