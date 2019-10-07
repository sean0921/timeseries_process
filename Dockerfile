FROM debian:buster
MAINTAINER holishing
WORKDIR /work
ADD . /work
RUN apt-get update \
 && apt-get upgrade -y \
 && apt-get install -y --no-install-recommends make automake autoconf autopoint autotools-dev gfortran gfortran-mingw-w64-i686 inkscape p7zip \
 && ./autogen.sh \
 && ./configure && make 7zpkg clean \
 && ./configure && make 7zpkg clean
