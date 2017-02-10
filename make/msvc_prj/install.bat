@echo off
rem set SRC=%~1%
rem set DST=%~2%
set SRC=..\..
set DST=..\..

echo ===========================================
echo Make directory tree

set INCDIR=%DST%\include
set SRCDIR=%DST%\src
rem set BINDIR=%DST%\bin
set SRCDIR_A=%DST%\altsrc
set DBGDIR_A=%DST%\dbglib
rem set LIBDIR_A=%DST%\lib

echo root folders
for %%d in ( %INCDIR% %SRCDIR% %BINDIR% %SRCDIR_A% %DBGDIR_A% %LIBDIR_A% ) do if not exist %%d (echo %%d & mkdir %%d)


set DIRLIST=access
set DIRLIST=%DIRLIST% algo
set DIRLIST=%DIRLIST% algo\blast
set DIRLIST=%DIRLIST% algo\blast\api
set DIRLIST=%DIRLIST% algo\blast\core
set DIRLIST=%DIRLIST% api
set DIRLIST=%DIRLIST% asnlib
set DIRLIST=%DIRLIST% biostruc
set DIRLIST=%DIRLIST% biostruc\cdd
set DIRLIST=%DIRLIST% biostruc\cn3d
set DIRLIST=%DIRLIST% cdromlib
set DIRLIST=%DIRLIST% cn3d
set DIRLIST=%DIRLIST% connect
set DIRLIST=%DIRLIST% connect\test
set DIRLIST=%DIRLIST% corelib
set DIRLIST=%DIRLIST% ctools
set DIRLIST=%DIRLIST% ddv
set DIRLIST=%DIRLIST% demo
set DIRLIST=%DIRLIST% desktop
set DIRLIST=%DIRLIST% gif
set DIRLIST=%DIRLIST% network\blast3\client
set DIRLIST=%DIRLIST% network\entrez\client
set DIRLIST=%DIRLIST% network\id1arch
set DIRLIST=%DIRLIST% network\medarch\client
set DIRLIST=%DIRLIST% network\nsclilib
set DIRLIST=%DIRLIST% network\spell\client
set DIRLIST=%DIRLIST% network\taxon1\common
set DIRLIST=%DIRLIST% network\taxon1\taxon2
set DIRLIST=%DIRLIST% network\vibnet
set DIRLIST=%DIRLIST% object
set DIRLIST=%DIRLIST% sequin
set DIRLIST=%DIRLIST% tools
set DIRLIST=%DIRLIST% util\creaders
set DIRLIST=%DIRLIST% util\tables
set DIRLIST=%DIRLIST% vibrant

echo src
for %%d in ( %DIRLIST% ) do if not exist %SRCDIR%\%%d (echo %SRCDIR%\%%d & mkdir %SRCDIR%\%%d)

echo ===========================================
echo Copy files

for %%d in ( %DIRLIST% ) do (echo %%d & if exist %SRC%\%%d copy %SRC%\%%d\*.h %INCDIR%)
for %%d in ( %DIRLIST% ) do (echo %%d & if exist %SRC%\%%d copy %SRC%\%%d\*.c %SRCDIR%\%%d)

for %%d in ( %DIRLIST% ) do (echo %%d & if exist %SRC%\%%d copy %SRC%\%%d\*.c %SRCDIR_A%)
for %%d in ( %DIRLIST% ) do (echo %%d & if exist %SRC%\%%d copy %SRC%\%%d\*.h %SRCDIR_A%)

set BUILDS=Debug
set BUILDS=%BUILDS% Release
set BUILDS=%BUILDS% DebugDLL
set BUILDS=%BUILDS% ReleaseDLL

echo build folders
for %%b in ( %BUILDS% ) do if not exist %DST%\%%b (echo %DST%\%%b & mkdir %DST%\%%b)

echo ===========================================
echo Copy files
for %%b in (%BUILDS%) do (if exist %SRC%\make\msvc_prj\corelib\ncbimain\%%b\ncbimain.obj  copy %SRC%\make\msvc_prj\corelib\ncbimain\%%b\ncbimain.obj %DST%\%%b)
for %%b in (%BUILDS%) do (if exist %SRC%\make\msvc_prj\corelib\ncbi\%%b\ncbithr.obj  copy %SRC%\make\msvc_prj\corelib\ncbi\%%b\ncbithr.obj %DST%\%%b)


set DIRLIST=%DIRLIST% asnlib\asntool
set DIRLIST=%DIRLIST% cdromlib\ncbiacc
set DIRLIST=%DIRLIST% cdromlib\ncbicdr
set DIRLIST=%DIRLIST% cdromlib\ncbinacc
set DIRLIST=%DIRLIST% cn3d\ncbicn3d
set DIRLIST=%DIRLIST% cn3d\ncbicn3d_ogl
set DIRLIST=%DIRLIST% corelib\ncbi
set DIRLIST=%DIRLIST% corelib\ncbimain
set DIRLIST=%DIRLIST% ddv\ddv
set DIRLIST=%DIRLIST% regexp
set DIRLIST=%DIRLIST% ddv\ddvlib
set DIRLIST=%DIRLIST% network\blast3\netblast
set DIRLIST=%DIRLIST% vibrant\vibrant
set DIRLIST=%DIRLIST% vibrant\vibrant_ogl

for %%d in ( %DIRLIST% ) do (for %%b in (%BUILDS%) do (if exist %SRC%\make\msvc_prj\%%d\%%b copy %SRC%\make\msvc_prj\%%d\%%b\*.lib %DST%\%%b))
for %%d in ( %DIRLIST% ) do (for %%b in (%BUILDS%) do (if exist %SRC%\make\msvc_prj\%%d\%%b copy %SRC%\make\msvc_prj\%%d\%%b\*.exe %DST%\%%b))

