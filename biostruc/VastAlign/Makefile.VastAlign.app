#
# $Id: Makefile.VastAlign.app,v 1.1 2005/07/26 15:44:12 chenj Exp $
#
# Author:  Jie Chen
#
# Build cgi "vastalign.cgi"
#
# $Log: Makefile.VastAlign.app,v $
# Revision 1.1  2005/07/26 15:44:12  chenj
# using lib center on frosty
#
#
#

REQUIRES = C-Toolkit PLATFORM_DB

APP = vastalign.cgi

SRC = cVAapp cVAcmd cVAres

LIB = pubvast shdb cjlib \
        entrez2 entrez2cli \
	$(PLATFORM_DB_LIB) \
	xcgi xcompress xconnect xcser xser xutil xncbi \
	$(Z_lib) 

CPPFLAGS = $(ORIG_CPPFLAGS) \
	$(NCBI_C_INCLUDE) \
        $(PLATFORM_DB_SYMBOLS) \
        -I$(includedir)/../src/internal/structure/DBLibs/VastSrv \
        -I$(includedir)/../src/internal/structure/DBLibs/Common \
	-I$(includedir)/../src/internal/structure/chenj \
	-I./

LIBS = $(ORIG_LIBS) \
	$(NCBI_C_LIBPATH) \
        -lncbimmdb -lncbitool -lncbiobj -lvibgif -lncbitxc2 -lnetcli \
	$(NCBI_C_ncbi) \
        $(PLATFORM_DB_LIBS) \
        $(Z_LIBS) $(DL_LIBS) $(NETWORK_LIBS)

