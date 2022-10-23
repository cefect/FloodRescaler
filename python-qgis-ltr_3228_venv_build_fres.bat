REM build a venv (need to run as admin)
SET VNAME=fres

call C:\LS\09_REPOS\01_COMMON\Qall\pyqgis_config\Q\3.22.8\python-qgis-ltr_3228.bat

SET VDIR=%userprofile%/.venv/%QVER%/%VNAME%

REM setup the virtual environment

ECHO building venv in %VDIR%
python -m venv --without-pip --copies "%VDIR%"
cd %VDIR%

ECHO launching pyvenv.cfg. set --system-site-packages = true
Notepad pyvenv.cfg
call "%VDIR%/Scripts/activate.bat"
cmd.exe

REM TODO: install from requirements?