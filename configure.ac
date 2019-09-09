#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([timeseries_process], [2019.909.0], [https://github.com/sean0921/timeseries_process])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
#AC_CONFIG_SRCDIR([])
#AC_CONFIG_HEADERS([config.h])
AC_PREFIX_DEFAULT(["$PWD"/build])
AC_PREFIX_PROGRAM(["$PWD"/build])

# Checks for programs.
AC_LANG(Fortran)
AC_PROG_FC
AC_PROG_F77
FCFLAGS="-std=legacy -g -static -cpp"
FFLAGS="-std=legacy -g -static -cpp"
AC_PROG_INSTALL

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([00_bern2time/Makefile
                 01_raw/Makefile
                 02_remove_trend/Makefile
                 03_remove_period/Makefile
                 04_remove_antenna/Makefile
                 05_Fit_velocity/Makefile
                 util/Makefile
                 Makefile])
AC_OUTPUT
