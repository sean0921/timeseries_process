## GNU makefile

TARGETdir=../build/$$(basename "$$PWD")
TARGET_PROGRAMS=remove_period
remove_period_SOURCES=remove_period.for

if IS_WINDOWS
remove_period_DEPENDENCIES=remove_period.rc.o
remove_period_LDFLAGS=remove_period.rc.o
SUFFIXES = .rc .o
RCFILES=remove_period.rc
RCOBJECT=remove_period.rc.o

$(RCFILES): %.rc:
	printf "id ICON \"$$(basename $@ .rc).ico\"" > $@
	cat ../common/version.rc >> $@

$(RCOBJECT): %.rc.o: %.rc
	convert -density 384 ../icons/hatch_wave.svg -define icon:auto-resize $$(basename $< .rc).ico
	$(WINDRES) $< -O coff -o $@

clean-local:
	rm -f *.ico *.rc

endif
