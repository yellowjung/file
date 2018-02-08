echo off
pushd "%~dp0"
cls

ver
ver > log.txt

echo.
echo Starting scan on C:\
echo.

find C:\* -type f -perm -u=r | dos2unix | file -f -