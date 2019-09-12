#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([timeseries_process], [2019.911.3], [https://github.com/sean0921/timeseries_process])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
#AC_CONFIG_SRCDIR([])
#AC_CONFIG_HEADERS([config.h])
AC_PREFIX_DEFAULT(["$PWD"/build])
AC_PREFIX_PROGRAM(["$PWD"/build])

# Checks for programs.
AC_LANG(Fortran)
AC_PROG_FC
AC_PROG_F77
FCFLAGS=" -std=legacy -g -static -cpp"
FFLAGS=" -std=legacy -g -static -cpp"
AC_CANONICAL_HOST
AM_PROG_AS
AC_PROG_INSTALL

## https://stackoverflow.com/questions/15398632
PLATF=""
AS_CASE([$host_os], [*mingw*], [PLATF="WINDOWS"])
AM_CONDITIONAL([IS_WINDOWS], [test "$PLATF" = "WINDOWS"])

if [[ $host = 'i686-w64-mingw32' ]] || [[ $host = 'x86_64-w64-mingw32' ]] || [[ $host = 'i686-w64-mingw32.static' ]] || [[ $host = 'x86_64-w64-mingw32.static' ]] || [[ $host = 'i686-w64-mingw32.shared' ]] || [[ $host = 'x86_64-w64-mingw32.shared' ]]
then
AC_SUBST(WINDRES,["$host_alias"-windres])
else
AC_SUBST(WINDRES,[windres])
fi

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
# Populate $host_cpu, $host_os, etc.
case $host_os in
  *mingw*)
FCFLAGS+=" -DMINGW"
FFLAGS+=" -DMINGW"
## https://stackoverflow.com/questions/7490978/
AC_CHECK_PROG(INKSCAPE_CHECK,inkscape,yes)
AS_IF([test x"$INKSCAPE_CHECK" != x"yes"], [AC_MSG_ERROR([Please install Inkscape before configuring.])])
AC_CHECK_PROG(IMAGEMAGICK_CHECK,convert,yes)
AS_IF([test x"$IMAGEMAGICK_CHECK" != x"yes"], [AC_MSG_ERROR([Please install imagemagick/graphicmagick before configuring.])])
    ;;
  *)
    ;;
esac

# Checks for library functions.

AC_CONFIG_FILES([00_bern2time/Makefile
                 01_raw/Makefile
                 02_remove_trend/Makefile
                 03_remove_period/Makefile
                 04_remove_antenna/Makefile
                 05_Fit_velocity/Makefile
                 util/Makefile
                 common/version.rc
                 Makefile])
AC_OUTPUT
