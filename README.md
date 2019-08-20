# Time Series Processing Tools

![](https://i.imgur.com/SmBl5UK.png)

## Compiling by GNU Fortran
* Linux / Other Unix-like platform:
    - Install `make`(only tested in GNU make), `gfortran`(maybe included in `gcc` in some OS/distribution).
    - type `make all` in each subfolder, it should be compiled (with debug symbol by default).

* Win32 platform
    - not test to built directly on Windows platform/MSYS2.
    - it is suggest to install `mingw-w64` cross-compile toolchain, for e.g., [`mxe`](https://mxe.cc) is a good choice.
    - mingw example in "Makefile"s in this repo is designed for [MXE](https://mxe.cc), you can adapted your mingw compiler command for that can fit your need.
    - type `make mingw` manully in each subfolder, it should be compiled (with debug symbol by default).
    - mingw build mode is NOT automatically detected, you should type `make mingw` on your own, so that the content of the program will be correct

# Temporary Install path
* Not really set yet, just make a temporary test directory named `test`, and place on there.
