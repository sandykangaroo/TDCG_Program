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

:Build
if exist "%BIN_ROOT%ipsxe-comp-vars.bat" @call "%BIN_ROOT%ipsxe-comp-vars.bat" intel64 vs2019

@echo on

ifort ./Code/*.f90 ./lib/*.lib ^
      /Od /exe:main.exe /assume:bscc /F0x1000000000 ^
      /module:./x64/Release/ /Fo:./x64/Release/
