# Time Series Processing Tools

![](https://i.imgur.com/JSbZLaM.png)

## Compiling by GNU Fortran and GNU Autotools

* Linux / Other Unix-like platform:
    - Install `make`(only tested in GNU make), `gfortran`(maybe included in `gcc` in some OS/distribution).
        + If you are using git version, you should install GNU Autotools(autotools-dev, aclocal, automake, autoconf), and remember to execute `./autogen.sh` first.
    - Type `./configure && make && make install` in main folder, it should be compiled (with debug symbol by default), and you can find the excutable on `build/` folder.

* Win32 platform
    - cross compiling by `mingw-w64` toolchain on Unix-like Platform (Linux, FreeBSD...etc) (**Faster**)
      + [MXE](https://mxe.cc) is a good choice, or you can search for `mingw-w64-*` prefix package by your package manager.
      + MinGW example in "Makefile"s in this repo is designed for [MXE](https://mxe.cc), you can adapted your MinGW compiler command by set `FC_MINGW` variable, for that can fit your need.
      + type cross-compile configure options like `./configure --build i686-pc-linux-gnu --host i686-w64-mingw32.static && make && make install` in main folder, it should be compiled (with debug symbol by default), and you can find the excutable on `build/` folder.
      + configure option format: `./configure <build platform prefix> <runtime platform prefix>`
      + If you are using git version, you should install GNU Autotools(autotools-dev, aclocal, automake, autoconf), and remember to execute `./autogen.sh` first.
    - compiling in Windows 7/10
      + [MSYS2](https://www.msys2.org) is a good choice, other tools are not tested or supported.
      + Type `./configure && make && make install` in main folder, it should be compiled (with debug symbol by default), and you can find the excutable on `build/` folder.
      + If you are using git version, you should install GNU Autotools(autotools-dev, aclocal, automake, autoconf), and remember to execute `./autogen.sh` first.

## Make a binary package (for easier to create tarballs)

* `make 7zpkg` (static binaries)
    - Please install `p7zip` package first.

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

## References/Related Articles

* Wdowinski, S., Bock, Y., Zhang, J., Fang, P., and Genrich, J. ( 1997), Southern California permanent GPS geodetic array: Spatial filtering of daily positions for estimating coseismic and postseismic displacements induced by the 1992 Landers earthquake, J. Geophys. Res., 102( B8), 18057– 18070, doi:10.1029/97JB01378.
* Nikolaidis, Rosanne. (2002). Observation of Geodetic and Seismic Deformation with the Global Positioning System.

