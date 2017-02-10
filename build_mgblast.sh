#!/bin/env bash
#customized build file to prevent building the GUI programs
os=$(uname -s) #Linux, Darwin
mtype=$(uname -m) # x86_64 or i686, i386 etc.
export HAVE_MOTIF=0
export HAVE_OGL=0
export HAVE_MAC=0

case "$os" in
  *Linux*)
     if [ "$mtype" == *x86_64* ]; then
       platform=linux64
     else
       platform=linux
     fi
     ;;
   FreeBSD)
     platform=freebsd
     ;;
   Darwin)
     platform=darwin
     ;;
   *)
     echo "Error: unrecognized platform $os" >&2
     exit 1
esac

echo "platform=$platform"
NCBI_DOT_MK="./platform/${platform}.ncbi.mk"
if [ ! -f "$NCBI_DOT_MK" ]; then
     echo "Error: cannot find platform file $NCBI_DOT_MK" >&2
     exit 1
fi

export NCBI_GNUTLS_INCLUDE=
export NCBI_GNUTLS_LIBS=

if  [ "$NCBI" == "1" ]; then
    if [ -f "$NCBI/ncbi.mk" ]; then
        unsplit='{ if (/\\$/) { l = l substr($0, 0, length-1) } else { print l $0; l="" } }'
        export NCBI_GNUTLS_INCLUDE=`awk "$unsplit" $NCBI/ncbi.mk | sed -ne 's/  */ /g; s/^ *NCBI_GNUTLS_INCLUDE *= *//p'`
        export NCBI_GNUTLS_LIBS=`awk "$unsplit" $NCBI/ncbi.mk | sed -ne 's/  */ /g; s/^ *NCBI_GNUTLS_LIBS *= *//p'`
    fi
else 
   if $(pkg-config gnutls --exists > /dev/null 2>&1) ; then
     export NCBI_GNUTLS_INCLUDE="`pkg-config gnutls --cflags` -DHAVE_LIBGNUTLS"
     export NCBI_GNUTLS_LIBS="`pkg-config gnutls --libs --static`"
   else 
     if $(libgnutls-config --version > /dev/null 2>&1) ; then
       export NCBI_GNUTLS_INCLUDE="`libgnutls-config --cflags` -DHAVE_LIBGNUTLS"
       export NCBI_GNUTLS_LIBS="`libgnutls-config --libs`"
     fi
   fi
fi

eval $(sed -e 's/^ *#.*//g' -e '/^$/d' -e 's/ *= */=/g' -e 's/^\([^=]*\)=\(.*\)$/export \1="\2"/' < $NCBI_DOT_MK)

export NCBI_OPTFLAG="-DNDEBUG $NCBI_OPTFLAG"

cd ./build
ln -s ../make/*.unx .
ln -s ../make/ln-if-absent .
mv makeall.unx makefile

if [ -n "$LIBPNG_DIR" ] &&  [ -n "$ZLIB_DIR" ]; then
    export PNG_INCLUDE="-D_PNG -I$LIBPNG_DIR -I$ZLIB_DIR"
    export PNG_LIBS="$LIBPNG_DIR/libpng.a $ZLIB_DIR/libz.a"
else
    export PNG_INCLUDE=""
    export PNG_LIBS=""
fi
export VIBWWWBLAST= 
export NONVIBWWWBLAST= 
export WWWBLAST=($VIBWWWBLAST $NONVIBWWWBLAST)
export OGL_NCBI_LIBS=""
export OGL_INCLUDE=""
export OGL_LIBS=""
export OGL_TARGETS=""

export ALL_VIB=
export DEMO_VIB=
export VIB="blastcl3 idfetch asn2gb tbl2asn gene2xml "
export NET_VIB=$VIB

export MFLG=""

CMD='make $MFLG GNUTLS_INCLUDE=\"$NCBI_GNUTLS_INCLUDE\" \
   GNUTLS_LIBS=\"$NCBI_GNUTLS_LIBS\" \
   CFLAGS1=\"$NCBI_OPTFLAG $NCBI_CFLAGS1 $OGL_INCLUDE $PNG_INCLUDE\" \
   LDFLAGS1=\"$NCBI_LDFLAGS1\" OTHERLIBS=\"$NCBI_OTHERLIBS\" \
   SHELL=\"$NCBI_MAKE_SHELL\" LCL=\"$NCBI_DEFAULT_LCL\" \
   RAN=\"$NCBI_RANLIB\" AR=\"$NCBI_AR\" CC=\"$NCBI_CC\" $ALL_VIB all'
eval echo $CMD
eval echo $CMD | sh 
make_stat=$?

if [ $make_stat -ne 0 ]; then
	cat <<EoF

Fatal error building NCBI core libraries.

EoF
	exit 1
fi

CMD='make $MFLG -f makeprog.unx GNUTLS_INCLUDE=\"$NCBI_GNUTLS_INCLUDE\" \
   GNUTLS_LIBS=\"$NCBI_GNUTLS_LIBS\" CFLAGS1=\"$NCBI_OPTFLAG $NCBI_CFLAGS1\" \
   LDFLAGS1=\"$NCBI_LDFLAGS1\" SHELL=\"$NCBI_MAKE_SHELL\" OTHERLIBS=\"$NCBI_OTHERLIBS\" \
   LCL=\"$NCBI_DEFAULT_LCL\" RAN=\"$NCBI_RANLIB\" AR=\"$NCBI_AR\" CC=\"$NCBI_CC\" $DEMO_VIB'
eval echo $CMD
eval echo $CMD | sh 

demo_stat=$?

# building other libs with multi-thread support
if [ -n "$MAKE_NET" ]; then
 if [ -n "$NCBI_MT_OTHERLIBS" ]; then
  CMD='make $MFLG -f makenet.unx GNUTLS_INCLUDE=\"$NCBI_GNUTLS_INCLUDE\" \
         GNUTLS_LIBS=\"$NCBI_GNUTLS_LIBS\" \
		 CFLAGS1=\"$NCBI_OPTFLAG $NCBI_CFLAGS1 $OGL_INCLUDE\" \
		 LDFLAGS1=\"$NCBI_LDFLAGS1\" SHELL=\"$NCBI_MAKE_SHELL\" \
		 AR=\"$NCBI_AR\" CC=\"$NCBI_CC\" RAN=\"$NCBI_RANLIB\" OTHERLIBS=\"$NCBI_OTHERLIBS\" \
		 THREAD_OBJ=$NCBI_THREAD_OBJ \
		 THREAD_OTHERLIBS=\"$NCBI_MT_OTHERLIBS\" \
		 NETENTREZVERSION=\"$NETENTREZVERSION\" $NET_VIB'
 else
  CMD='make $MFLG -f makenet.unx GNUTLS_INCLUDE=\"$NCBI_GNUTLS_INCLUDE\" \
         GNUTLS_LIBS=\"$NCBI_GNUTLS_LIBS\" \
		 CFLAGS1=\"$NCBI_OPTFLAG $NCBI_CFLAGS1 $OGL_INCLUDE\" \
		 LDFLAGS1=\"$NCBI_LDFLAGS1\" SHELL=\"$NCBI_MAKE_SHELL\" \
		 AR=\"$NCBI_AR\" CC=\"$NCBI_CC\" RAN=\"$NCBI_RANLIB\" OTHERLIBS=\"$NCBI_OTHERLIBS\" \
		 NETENTREZVERSION=\"$NETENTREZVERSION\" $NET_VIB'
 fi
 eval echo $CMD
 eval echo $CMD | sh 

fi

net_stat=$?

if [ "$demo_stat" -ne 0 ] || [ "$net_stat" -ne 0 ]; then
 echo "Error encountered during build!"
 exit 1
fi

echo '*********************************************************'
echo '*The new binaries are located in ./build/ directory*'
echo '*********************************************************'
