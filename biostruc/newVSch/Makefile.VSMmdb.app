#
# $Id: Makefile.VSMmdb.app,v 1.1 2005/07/26 17:05:42 chenj Exp $
#
# Author: Jie Chen
# Build cgi for new VS
#
# $Log: Makefile.VSMmdb.app,v $
# Revision 1.1  2005/07/26 17:05:42  chenj
# Making linux version cgi
#
#
#
#
#


REQUIRES = C-Toolkit PLATFORM_DB

APP = VSMmdb2.cgi

SRC = cVSMmdbapp cVSMmdbcmd cVSMmdbres mmdbuti SendSummary qmuti

#LIB = pubstruct qman dart shdb cjlib \

#LIB = cjlib \

LIB = qman shdb cjlib \
  	seq pub medline biblio general taxon1 \
        entrez2 entrez2cli \
	$(PLATFORM_DB_LIB) \
	xcgi xhtml xcompress xconnect xcser xser xutil xncbi \
	$(Z_lib) 

VSLibDir = /web/public/htdocs/Structure/chenj/vast/VSchDB/DBLibs
HeDir = /netopt/structure/lib
PRE_LIBS = -L$(VSLibDir) -lVastSrch

	# -L$(HeDir) -lpubstruct -lqman -ldart -lshdb



CPPFLAGS = $(ORIG_CPPFLAGS) \
	$(NCBI_C_INCLUDE) \
        $(PLATFORM_DB_SYMBOLS) \
        -I$(includedir)/../src/internal/structure/DBLibs/VastSrv \
        -I$(includedir)/../src/internal/structure/DBLibs/MmdbSrv \
        -I$(includedir)/../src/internal/structure/DBLibs/Common \
	-I$(includedir)/../src/internal/structure/chenj \
	-I$(includedir)/../src/internal/structure/QueueManager \
	-I$(includedir)/../src/internal/structure/Dart/LIB \
	-I$(VSLibDir) \
	-I./

LIBS = $(ORIG_LIBS) \
	$(NCBI_C_LIBPATH) \
 	-lncbimmdb -lncbitool -lncbiobj -lvibgif -lncbitxc2 -lnetcli \
	$(NCBI_C_ncbi) \
        $(PLATFORM_DB_LIBS) \
        $(Z_LIBS) $(DL_LIBS) $(NETWORK_LIBS)

