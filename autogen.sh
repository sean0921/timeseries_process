#!/bin/sh
set -eux

test -e $(which aclocal)  || echo "Please install aclocal  !!"
test -e $(which automake) || echo "Please install automake !!"
test -e $(which autoconf) || echo "Please install autoconf !!"

aclocal \
&& automake --add-missing \
&& autoconf
