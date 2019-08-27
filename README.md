# Time Series Processing Tools

![](https://i.imgur.com/SmBl5UK.png)

## Compiling by GNU Fortran

* Linux / Other Unix-like platform:
    - Install `make`(only tested in GNU make), `gfortran`(maybe included in `gcc` in some OS/distribution).
    - Type `make all` in each subfolder, it should be compiled (with debug symbol by default).

* Win32 platform
    - Not tested yet to built directly on Windows platform/MSYS2.
    - It is suggest to install `mingw-w64`(for i686 or x86\_64) cross-compile toolchain, for e.g., [`mxe`](https://mxe.cc) is a good choice.
    - MinGW example in "Makefile"s in this repo is designed for [MXE](https://mxe.cc), you can adapted your MinGW compiler command for that can fit your need.
    - type `make mingw` manully in each subfolder, it should be compiled (with debug symbol by default).
    - MinGW build mode is NOT automatically detected, you should type `make mingw` on your own, so that the content of the program will be correct

## Input file

* `bern2time`:
    - `getpl.gout` (from `getpl`)
    - `sta-file` (from `getpl`)

* `comfilt` (`raw`
    - `comfilt.inp` (`raw.inp`, adapted from `sta-file`, format example:`comfilt.inp.example`/`raw.inp.example`)

* `remove_{trend,period}`
    - `remove.inp` (copied from `comfilt.inp`/`raw.inp`)

* `remove_attenna`:
    - `remove.inp` (copied from `comfilt.inp`/`raw.inp`)

* `fit`:
    - `fit.inp` (copied from `comfilt.inp`/`raw.inp`)

## make a binary package (tarball)

* `make package` (static linux binary)
* `make mingw_package` (static win32 excutable)
