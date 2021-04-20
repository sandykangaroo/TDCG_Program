@echo off

@REM Changed by Xueliang Li in 2021.4.18.
@REM Change file root in line 12 if computer changed.
@REM Root name can be find in:
@REM Win Start
@REM Intel Parallel Studio XE
@REM Shortcut file name "Compiler 19.0 Update 5 for Intel 64 Visual Studio 2019 environment"
@REM Right click
@REM Properties

:Start

:: Cache some environment variables.
set BIN_ROOT=C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2019.5.281\windows\bin\

set TARGET_ARCH=intel64

:Build

if defined INSPECTOR_2019_DIR (
    if exist "%INSPECTOR_2019_DIR%\inspxe-vars.bat" @call "%INSPECTOR_2019_DIR%\inspxe-vars.bat" quiet
)
if defined VTUNE_AMPLIFIER_2019_DIR (
    if exist "%VTUNE_AMPLIFIER_2019_DIR%\amplxe-vars.bat" @call "%VTUNE_AMPLIFIER_2019_DIR%\amplxe-vars.bat" quiet
)
if defined ADVISOR_2019_DIR (
    if exist "%ADVISOR_2019_DIR%\advixe-vars.bat" @call "%ADVISOR_2019_DIR%\advixe-vars.bat" quiet
)

if /i "%TARGET_ARCH%"=="intel64" (
:: ITAC
    if exist "%BIN_ROOT%..\..\..\itac_2019\bin\itacvars.bat" (
        @call "%BIN_ROOT%..\..\..\itac_2019\bin\itacvars.bat"
    ) else (
        if exist "%BIN_ROOT%..\..\..\..\itac_2019\bin\itacvars.bat" @call "%BIN_ROOT%..\..\..\..\itac_2019\bin\itacvars.bat"
    )
)

if exist "%BIN_ROOT%compilervars.bat" @call "%BIN_ROOT%compilervars.bat" intel64

@echo on

ifort ./Code/*.f90 ./lib/*.lib ^
      -O3 /exe:main.exe /assume:bscc /F0x1000000000 ^
      /module:./x64/Release/ /Fo:./x64/Release/
