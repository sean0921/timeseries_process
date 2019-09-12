## GNU makefile

TARGETdir=../build/$$(basename "$$PWD")
TARGET_PROGRAMS=comfilt raw comfilt_trend
comfilt_SOURCES=comfilt.for ../common/bubblesort.for ../common/int2char.for ../common/lsqr.for
raw_SOURCES=comfilt.for ../common/bubblesort.for ../common/int2char.for ../common/lsqr.for
raw_FCFLAGS=-DRAW
raw_FFLAGS=-DRAW
comfilt_trend_SOURCES=comfilt_trend.for ../common/bubblesort.for ../common/int2char.for ../common/lsqr.for
TARGET_DATA=comfilt.inp.example  raw.inp.example

if IS_WINDOWS
comfilt_DEPENDENCIES=comfilt.rc.o
raw_DEPENDENCIES=raw.rc.o
comfilt_trend_DEPENDENCIES=comfilt_trend.rc.o
comfilt_LDFLAGS=comfilt.rc.o
raw_LDFLAGS=raw.rc.o
comfilt_trend_LDFLAGS=comfilt_trend.rc.o
SUFFIXES = .rc .o
RCFILES=comfilt.rc raw.rc comfilt_trend.rc
RCOBJECT=comfilt.rc.o comfilt_trend.rc.o

$(RCFILES): %.rc:
	printf "id ICON \"$$(basename $@ .rc).ico\"" > $@
	cat ../common/version.rc >> $@

comfilt.rc.o: %.rc.o: %.rc
	convert -density 384 ../icons/lines-and-angles_filter.svg -define icon:auto-resize $$(basename $< .rc).ico
	$(WINDRES) $< -O coff -o $@

raw.rc.o: %.rc.o: %.rc
	convert -density 384 ../icons/lines-and-angles_adjust.svg -define icon:auto-resize $$(basename $< .rc).ico
	$(WINDRES) $< -O coff -o $@

comfilt_trend.rc.o: %.rc.o: %.rc
	convert -density 384 ../icons/lines-and-angles_chart.svg -define icon:auto-resize $$(basename $< .rc).ico
	$(WINDRES) $< -O coff -o $@

clean-local:
	rm -f *.ico *.rc

endif
