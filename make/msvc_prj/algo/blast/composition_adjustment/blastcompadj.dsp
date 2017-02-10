# Microsoft Developer Studio Project File - Name="blastcompadj" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=blastcompadj - Win32 DebugDLL
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "blastcompadj.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "blastcompadj.mak" CFG="blastcompadj - Win32 DebugDLL"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "blastcompadj - Win32 DebugDLL" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe
# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "DebugDLL"
# PROP BASE Intermediate_Dir "DebugDLL"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "DebugDLL"
# PROP Intermediate_Dir "DebugDLL"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /I "..\..\.." /I "..\..\..\corelib" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /I "..\..\..\..\.." /I "..\..\..\..\..\corelib" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
# Begin Target

# Name "blastcompadj - Win32 DebugDLL"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\compo_heap.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\compo_mode_condition.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\composition_adjustment.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\matrix_frequency_data.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\nlm_linear_algebra.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\optimize_target_freq.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\redo_alignment.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\smith_waterman.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\compo_heap.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\compo_mode_condition.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\composition_adjustment.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\composition_constants.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\matrix_frequency_data.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\nlm_linear_algebra.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\optimize_target_freq.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\redo_alignment.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\algo\blast\composition_adjustment\smith_waterman.h
# End Source File
# End Group
# End Target
# End Project
