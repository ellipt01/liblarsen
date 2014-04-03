#!/bin/sh

if [ -e Makefile ]; then
	echo "make clean"
	make clean > /dev/null 2>&1
fi

echo "rm -rf autom4te.cache"
rm -rf autom4te.cache

echo "rm -f config.h.in install-sh missing depcomp aclocal.m4 acinclude.m4"
rm -f config.h.in install-sh missing depcomp aclocal.m4 acinclude.m4

echo "rm -f stamp-h.in"
rm -f stamp-h.in

echo "rm -f ltmain.sh stamp-h1 config.log config.h config.sub config.status config.guess"
rm -f ltmain.sh stamp-h1 config.log config.h config.sub config.status config.guess

echo "rm -f Makefile include/Makefile src/Makefile"
rm -f Makefile include/Makefile src/Makefile

echo "rm -f configure libtool"
rm -f configure libtool

echo "rm -f include/Makefile.am.bak"
rm -f include/Makefile.am.bak

echo "rm -f src/Makefile.am.bak"
rm -f src/Makefile.am.bak

echo "rm -f Makefile.in src/Makefile.in include/Makefile.in"
rm -f Makefile.in src/Makefile.in include/Makefile.in

echo "rm -f src/.deps/*.Plo"
rm -f src/.deps/*.Plo

echo -n "rm "
find . -name '*~' -print
find . -name '*~' -exec rm {} \;

DEBUG=Debug
PAR=$(readlink -f $DEBUG)
if [ -d $PAR ]; then
	echo "cd $PAR"
	cd $PAR
	echo "make clean"
	make clean > /dev/null 2>&1
fi
