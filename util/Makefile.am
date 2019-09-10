## GNU makefile

TARGETdir=../build/$$(basename "$$PWD")
TARGET_PROGRAMS=plot_ts
plot_ts_SOURCES=plot_ts.for
TARGET_DATA=plot_ts.inp.example

if IS_WINDOWS
plot_ts_DEPENDENCIES=plot_ts.rc.o
plot_ts_LDFLAGS=plot_ts.rc.o
SUFFIXES = .rc .o
RCFILES=plot_ts.rc
RCOBJECT=plot_ts.rc.o

$(RCFILES): %.rc:
	printf "id ICON \"$$(basename $@ .rc).ico\"" > $@
	cat ../common/version.rc >> $@

$(RCOBJECT): %.rc.o: %.rc
	convert -density 384 ../icons/sketchy_plot.svg -define icon:auto-resize $$(basename $< .rc).ico
	$(WINDRES) $< -O coff -o $@

clean-local:
	rm -f *.ico *.rc

endif
