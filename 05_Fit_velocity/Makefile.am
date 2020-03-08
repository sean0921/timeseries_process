## GNU makefile

TARGETdir=../build/$$(basename "$$PWD")
TARGET_PROGRAMS=fit
fit_SOURCES=fit.for ../common/int2char.for ../common/lsqr.for

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
	$(WINDRES) $< -O coff -o $@

ico: %.rc
	convert -density 384 ../icons/dashing_plot.svg -define icon:auto-resize $$(basename $< .rc).ico

clean-local:
	rm -f *.rc

endif
