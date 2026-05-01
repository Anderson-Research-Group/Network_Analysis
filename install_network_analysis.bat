REM remove existing network_analysis env if any
set CONDAENV=network_analysis
IF NOT "%1"=="" SET CONDAENV=%1
echo "creating new conda environment for Network Analysis called %CONDAENV%..."

REM remove any existing environment with this name
call conda env remove -n %CONDAENV%

REM update anaconda
call conda config --add channels anaconda
call conda config --add channels conda-forge

REM install network_analysis deps
call conda create --yes --name %CONDAENV% python=3.9 pip=25.2 openssl==3.6.0 || goto :error
call conda activate %CONDAENV%

REM different versions of open3d for different OSes, so we install it manually here
call pip install open3d==0.19.0 || goto :error

REM reactivate network_analysis environment
call conda activate base
call conda activate %CONDAENV%

call python -m pip install -r python_requirements.txt || goto :error

echo "Network Analysis environment %CONDAENV% is ready to use."
goto :EOF

:error
echo "Installation failed. Please check the error message above."
exit /b %ERRORLEVEL%