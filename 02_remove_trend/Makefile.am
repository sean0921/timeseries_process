## GNU makefile

TARGETdir=../build/$$(basename "$$PWD")
TARGET_PROGRAMS=remove_trend
remove_trend_SOURCES=remove_trend.for

if IS_WINDOWS
remove_trend_DEPENDENCIES=remove_trend.rc.o
remove_trend_LDFLAGS=remove_trend.rc.o
SUFFIXES = .rc .o
RCFILES=remove_trend.rc
RCOBJECT=remove_trend.rc.o

$(RCFILES): %.rc:
	printf "id ICON \"$$(basename $@ .rc).ico\"" > $@
	cat ../common/version.rc >> $@

$(RCOBJECT): %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

ico: %.rc
	convert -density 384 ../icons/lines-and-angles_trend.svg -define icon:auto-resize $$(basename $< .rc).ico

clean-local:
	rm -f *.rc

endif
