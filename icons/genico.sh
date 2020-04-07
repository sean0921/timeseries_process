#!/bin/sh
set -eux

for i in $(ls *.svg)
do
    convert -density 384 $i -define icon:auto-resize $(basename $i .svg).ico
done
