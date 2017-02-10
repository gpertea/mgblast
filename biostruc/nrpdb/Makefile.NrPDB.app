#
# $Id: Makefile.NrPDB.app,v 1.1 2005/07/26 17:27:00 chenj Exp $
#
# Author:  Jie Chen
# Build nrpdbsrv.cgi
#
# $Log#
#
#

REQUIRES = C-Toolkit PLATFORM_DB

APP = nrpdbsrv.cgi

SRC = nrpdbsrv

LIB = xcgi xhtml xcompress xconnect xcser xser xutil xncbi \
	$(Z_lib) 

CPPFLAGS = $(ORIG_CPPFLAGS) $(NCBI_C_INCLUDE)

LIBS = $(ORIG_LIBS) \
	$(NCBI_C_LIBPATH) \
        $(NCBI_C_ncbi)

