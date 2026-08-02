@echo off
setlocal EnableExtensions EnableDelayedExpansion

:: Always work in the folder containing this batch file.
cd /d "%~dp0"

:: Use the existing Anaconda Python installation.
set "PYTHON=D:\Anaconda3\python.exe"
set "CONVERTER=%~dp0html_to_pdf.py"

if not exist "%PYTHON%" (
    echo Python not found: "%PYTHON%"
    exit /b 1
)

if not exist "%CONVERTER%" (
    echo Converter not found: "%CONVERTER%"
    exit /b 1
)

:: Search for an HTML file whose name contains 0.h by default.
set "defaultPattern=0.h"
set "pattern=!defaultPattern!"

:: Select the first matching HTML file, in the same order as "convert to html.bat".
set "found="
for /f "delims=" %%i in ('dir /b /a-d "*!pattern!*.html" 2^>nul') do (
    set "found=%%i"
    goto :default_found
)

:default_found
if not defined found (
    echo No HTML file matching "!defaultPattern!" was found.
    goto :manual_input
)

echo Default match: "!found!"
choice /c YN /n /m "Use this HTML file? [Y/N] "
if errorlevel 2 goto :manual_input
if errorlevel 1 goto :run

:manual_input
set "pattern="
set /p "pattern=Enter a keyword to match an HTML file: "
if not defined pattern (
    echo No keyword entered.
    exit /b 1
)

set "found="
for /f "delims=" %%i in ('dir /b /a-d "*!pattern!*.html" 2^>nul') do (
    set "found=%%i"
    goto :manual_found
)

:manual_found
if not defined found (
    echo No HTML file matching "!pattern!" was found.
    exit /b 1
)

echo Found: "!found!"
choice /c YN /n /m "Use this HTML file? [Y/N] "
if errorlevel 2 (
    echo Cancelled.
    exit /b 1
)

:run
"%PYTHON%" "%CONVERTER%" "%found%"
echo.
echo Press any key to exit...
pause >nul

endlocal
