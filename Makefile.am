#### Simple Makefile for build a pack

SUBDIRS = src

MAINdir=build/
MAIN_DATA=COPYRIGHT.md README.md

all:
	for i in ${SUBDIRS};do cd $$i;${MAKE};cd -;done

clean-local:
	for i in ${SUBDIRS};do cd $$i;${MAKE} clean;cd -;done

7zpkg: install
	p7zip build
	(mv -f build.7z ${PACKAGE}-$$(grep Target config.log |awk '{print $$2}'|uniq)-${PACKAGE_VERSION}.7z) || mv -f build.7z ${PACKAGE}-${PACKAGE_VERSION}.7z

if IS_WINDOWS
MAIN_DATA+= LICENSE_icons.txt
windres:
	@echo @WINDRES@
endif
