# Microsoft Developer Studio Project File - Name="ncbitxc2" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=ncbitxc2 - Win32 DebugDLL
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "ncbitxc2.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ncbitxc2.mak" CFG="ncbitxc2 - Win32 DebugDLL"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ncbitxc2 - Win32 DebugDLL" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "ncbitxc2 - Win32 DebugDLL"

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
# ADD BASE CPP /nologo /MDd /W3 /Gm /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# SUBTRACT BASE CPP /Fr
# ADD CPP /nologo /MDd /W3 /Gm /GR /Z7 /Od /I "..\..\..\..\.." /I "..\..\..\..\..\asnlib" /I "..\..\..\..\..\api" /I "..\..\..\..\..\corelib" /I "..\..\..\..\..\object" /I "..\..\..\..\..\network\taxon1\common" /I "..\..\..\..\..\network\taxon1\taxon2" /I "..\..\..\..\..\network\nsclilib" /D "WIN32" /D "_DEBUG" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "ncbitxc2 - Win32 DebugDLL"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\..\..\network\taxon1\common\objtax1.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\network\taxon1\taxon2\parttree.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\network\taxon1\taxon2\taxcs.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\network\taxon1\common\taxutil.c
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\network\taxon1\taxon2\tc2proc.c
# ADD CPP /D "TAXSERVICE"
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\network\taxon1\taxon2\txcdproc.c
# ADD CPP /D "TAXSERVICE"
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\..\..\network\taxon1\common\objtax1.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\..\network\taxon1\common\taxutil.h
# End Source File
# End Group
# End Target
# End Project
