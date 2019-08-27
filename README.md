# Time Series Processing Tools

![](https://i.imgur.com/SmBl5UK.png)

## Compiling by GNU Fortran

* Linux / Other Unix-like platform:
    - Install `make`(only tested in GNU make), `gfortran`(maybe included in `gcc` in some OS/distribution).
    - Type `make all` in each subfolder, it should be compiled (with debug symbol by default).

* Win32 platform
    - cross compiling by `mingw-w64` toolchain on Unix-like Platform (Linux, FreeBSD...etc) (**Faster**)
      + [MXE](https://mxe.cc) is a good choice, or you can search for `mingw-w64-*` prefix package by your package manager.
      + MinGW example in "Makefile"s in this repo is designed for [MXE](https://mxe.cc), you can adapted your MinGW compiler command by set `FC_MINGW` variable, for that can fit your need.
      + type `make mingw` in main folder, or manully in each subfolder, it should be compiled (with debug symbol by default).
    - compiling in Windows 7/10
      + [MSYS2](https://www.msys2.org) is a good choice, other tools are not tested or supported.
      + You should adapted your MinGW compiler command by set `FC_MINGW` variable, for that can fit your need.
      + type `make mingw` in main folder, or manully in each subfolder, it should be compiled (with debug symbol by default).

## Input file

* `bern2time`:
    - `getpl.gout` (from `getpl`)
    - `sta-file` (from `getpl`)

* `comfilt` (`raw`)
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
