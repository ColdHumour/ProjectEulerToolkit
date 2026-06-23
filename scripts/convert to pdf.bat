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

set /p "pattern=Enter a keyword to match an HTML file: "
if "%pattern%"=="" (
    echo No keyword entered.
    exit /b 1
)

:: Select the first matching HTML file, in the same order as "convert to html.bat".
set "found="
for /f "delims=" %%i in ('dir /b /a-d "*%pattern%*.html" 2^>nul') do (
    set "found=%%i"
    goto :found_one
)

:found_one
if "%found%"=="" (
    echo No HTML file matching "%pattern%" was found.
    exit /b 1
)

echo Found: "%found%"
"%PYTHON%" "%CONVERTER%" "%found%"

endlocal
