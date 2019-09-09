#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([timeseries_process], [2019.909.0], [sean.li.shin.ho@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
#AC_CONFIG_SRCDIR([])
#AC_CONFIG_HEADERS([config.h])

AC_LANG(Fortran)
# Checks for programs.
AC_PROG_FC
FCFLAGS="-std=legacy -g -static -cpp"
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
