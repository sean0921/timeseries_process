## GNU makefile

TARGETdir=../build/$$(basename "$$PWD")
TARGET_PROGRAMS=fit
fit_SOURCES=fit.for

if IS_WINDOWS
fit_DEPENDENCIES=fit.rc.o
fit_LDFLAGS=fit.rc.o
SUFFIXES = .rc .o
RCFILES=fit.rc
RCOBJECT=fit.rc.o

$(RCFILES): %.rc:
	printf "id ICON \"$$(basename $@ .rc).ico\"" > $@
	cat ../common/version.rc >> $@

$(RCOBJECT): %.rc.o: %.rc
	convert -density 384 ../icons/dashing_plot.svg -define icon:auto-resize $$(basename $< .rc).ico
	$(WINDRES) $< -O coff -o $@

clean-local:
	rm -f *.ico *.rc

endif
