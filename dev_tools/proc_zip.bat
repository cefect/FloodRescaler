:: zipping processing directory
SET TARGET_ZIP=floodrescaler.zip

:: setup environment
call %~dp0../env/activate_py.bat

ECHO on

SET PATH=C:\ProgramData\chocolatey\bin;%PATH%

:: change to plugin directory
cd %SRC_DIR%
 
 

REM ** 2) Create a zip file **
 
if exist %TARGET_ZIP% del %TARGET_ZIP%
 
REM 7z a %TARGET_ZIP% %SRC_DIR%\floodrescaler\processing\
7z a %TARGET_ZIP% %SRC_DIR%\floodrescaler\processing\ -xr!/.pytest_cache -xr!/__pycache__ -xr!__init__.py 

 

REM ** 3) Copy the zip file **
::move cancurve.zip "%SRC_DIR%\plugin_zips" /Y
xcopy %TARGET_ZIP% %SRC_DIR%\deploy /Y /I

if exist %TARGET_ZIP% del %TARGET_ZIP%
 

cmd.exe /k