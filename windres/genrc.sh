#!/bin/sh
set -eux

## POSIX shell script for genenrate files for resources compilers (e.g. windres)

for i in bern2time comfilt comfilt_trend fit getpl plot_ts raw remove_antenna remove_period remove_trend
do cat > "$i".rc << EOF
#include "version.rc"
id ICON "$i.ico"
EOF
done
