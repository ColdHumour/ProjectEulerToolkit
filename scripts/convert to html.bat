@echo off
setlocal EnableExtensions EnableDelayedExpansion

:: Always search in the folder containing this batch file.
cd /d "%~dp0"

:: Change this command prefix if needed.
set "prefix=jupyter nbconvert --to html --template classic"

:: Search for a file whose name contains 0.i by default.
set "defaultPattern=0.i"
set "pattern=!defaultPattern!"

:: Search files in this folder; the first result is selected alphabetically.
set "found="
for /f "delims=" %%i in ('dir /b /a-d "*!pattern!*" 2^>nul') do (
    set "found=%%i"
    goto :default_found
)

:default_found
if not defined found (
    echo No file containing "!defaultPattern!" was found.
    goto :manual_input
)

echo Default match: "!found!"
choice /c YN /n /m "Use this file? [Y/N] "
if errorlevel 2 goto :manual_input
if errorlevel 1 goto :run

:manual_input
set "pattern="
set /p "pattern=Enter a keyword to match a file: "
if not defined pattern (
    echo No keyword entered.
    exit /b 1
)

set "found="
for /f "delims=" %%i in ('dir /b /a-d "*!pattern!*" 2^>nul') do (
    set "found=%%i"
    goto :manual_found
)

:manual_found
if not defined found (
    echo No file containing "!pattern!" was found.
    exit /b 1
)

echo Found: "!found!"
choice /c YN /n /m "Use this file? [Y/N] "
if errorlevel 2 (
    echo Cancelled.
    exit /b 1
)

:: Run the conversion command with the selected file.
:run
%prefix% "!found!"
echo.
echo Press any key to exit...
pause >nul

endlocal
