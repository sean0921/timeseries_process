## GNU makefile

TARGETdir=../build/$$(basename "$$PWD")
TARGET_PROGRAMS=bern2time getpl
bern2time_SOURCES=bern2time.for ../common/bubblesort.for
getpl_SOURCES=getpl.for ../common/bubblesort.for

if IS_WINDOWS
bern2time_DEPENDENCIES=bern2time.rc.o
getpl_DEPENDENCIES=getpl.rc.o
bern2time_LDFLAGS=bern2time.rc.o
getpl_LDFLAGS=getpl.rc.o
SUFFIXES = .rc .o
RCFILES=bern2time.rc getpl.rc

$(RCFILES): %.rc:
	printf "id ICON \"$$(basename $@ .rc).ico\"" > $@
	cat ../common/version.rc >> $@

bern2time.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

getpl.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

ico:
	convert -density 384 ../icons/lines-and-angles_report.svg -define icon:auto-resize $$(basename $< .rc).ico
	convert -density 384 ../icons/lines-and-angles_locate.svg -define icon:auto-resize $$(basename $< .rc).ico

clean-local:
	rm -f *.rc

endif
