#! /bin/sh

echo "If you use a recent version of autotools, this script is obsolete"
echo "Just run autoreconf -i -f -v"
echo "followed by ./configure"
echo 

# Run this to generate all the auto-generated files needed by the GNU
# configure program
libtoolize --automake
aclocal
autoheader
automake --add-missing --gnu --force-missing
autoconf
echo "Now use ./configure"
