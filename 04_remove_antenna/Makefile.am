## GNU makefile

TARGETdir=../build/$$(basename "$$PWD")
TARGET_PROGRAMS=remove_antenna
remove_antenna_SOURCES=remove_antenna.for

if IS_WINDOWS
remove_antenna_DEPENDENCIES=remove_antenna.rc.o
remove_antenna_LDFLAGS=remove_antenna.rc.o
SUFFIXES = .rc .o
RCFILES=remove_antenna.rc
RCOBJECT=remove_antenna.rc.o

$(RCFILES): %.rc:
	printf "id ICON \"$$(basename $@ .rc).ico\"" > $@
	cat ../common/version.rc >> $@

$(RCOBJECT): %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

ico: %.rc
	convert -density 384 ../icons/science-and-fiction_listen.svg -define icon:auto-resize $$(basename $< .rc).ico

clean-local:
	rm -f *.rc

endif
