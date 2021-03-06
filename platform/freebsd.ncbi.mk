#
# $Id: freebsd.ncbi.mk,v 1.9 2003/05/07 17:40:12 beloslyu Exp $
#
# This platform is not supported. Tested on FreeBSD 4.5-RELEASE
#
NCBI_DEFAULT_LCL = fbd
NCBI_MAKE_SHELL = /bin/sh
NCBI_AR=ar
NCBI_CC = gcc -pipe -pthread
NCBI_CFLAGS1 = -c
NCBI_LDFLAGS1 = -O
NCBI_OPTFLAG = -O
NCBI_BIN_MASTER = /home/coremake/ncbi/bin
NCBI_BIN_COPY = /home/coremake/ncbi/bin
NCBI_INCDIR = /home/coremake/ncbi/include
NCBI_LIBDIR = /home/coremake/ncbi/lib
NCBI_ALTLIB = /home/coremake/ncbi/altlib
#will work only when you have Motif installed!
NCBI_VIBFLAG = -I/usr/X11R6/include -L/usr/X11R6/lib -DWIN_MOTIF
NCBI_VIBLIBS = -lXm -lXmu -lXp -lXpm -lXt -lX11 -lXext 
NCBI_DISTVIBLIBS = -L/usr/X11R6/lib /usr/X11R6/lib/libXm.a /usr/X11R6/lib/libXp.a /usr/X11R6/lib/libXpm.a -lXmu -lXt -lX11 -lXext
NCBI_OTHERLIBS = -lm
NCBI_RANLIB = ranlib
# Used by makedis.csh
NCBI_MT_OTHERLIBS = -pthread
NCBI_OTHERLIBS_MT = $(NCBI_MT_OTHERLIBS) -lm
NCBI_THREAD_OBJ = ncbithr.o
NETENTREZVERSION = 2.02c2ASN1SPEC6 

NCBI_LBSM_SRC = ncbi_lbsmd_stub.c
NCBI_LBSM_OBJ = ncbi_lbsmd_stub.o

