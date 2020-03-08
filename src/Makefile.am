#Makefile

TARGETdir=../build
TARGET_PROGRAMS=bern2time getpl comfilt raw comfilt_trend remove_trend remove_period remove_antenna fit plot_ts
TARGET_DATA=comfilt.inp.example  raw.inp.example plot_ts.inp.example order.txt

bern2time_SOURCES=bern2time.for ../common/bubblesort.for
getpl_SOURCES=getpl.for ../common/bubblesort.for
comfilt_SOURCES=comfilt.for ../common/bubblesort.for ../common/int2char.for ../common/lsqr.for
raw_SOURCES=comfilt.for ../common/bubblesort.for ../common/int2char.for ../common/lsqr.for
raw_FCFLAGS=-DRAW
raw_FFLAGS=-DRAW
comfilt_trend_SOURCES=comfilt_trend.for ../common/bubblesort.for ../common/int2char.for ../common/lsqr.for
remove_trend_SOURCES=remove_trend.for
remove_period_SOURCES=remove_period.for
remove_antenna_SOURCES=remove_antenna.for
fit_SOURCES=fit.for ../common/int2char.for ../common/lsqr.for
plot_ts_SOURCES=plot_ts.for

if IS_WINDOWS

bern2time_DEPENDENCIES=bern2time.rc.o
getpl_DEPENDENCIES=getpl.rc.o
comfilt_DEPENDENCIES=comfilt.rc.o
raw_DEPENDENCIES=raw.rc.o
comfilt_trend_DEPENDENCIES=comfilt_trend.rc.o
remove_trend_DEPENDENCIES=remove_trend.rc.o
remove_period_DEPENDENCIES=remove_period.rc.o
remove_antenna_DEPENDENCIES=remove_antenna.rc.o
fit_DEPENDENCIES=fit.rc.o
plot_ts_DEPENDENCIES=plot_ts.rc.o
bern2time_LDFLAGS=bern2time.rc.o
getpl_LDFLAGS=getpl.rc.o
comfilt_LDFLAGS=comfilt.rc.o
raw_LDFLAGS=raw.rc.o
comfilt_trend_LDFLAGS=comfilt_trend.rc.o
remove_trend_LDFLAGS=remove_trend.rc.o
remove_period_LDFLAGS=remove_period.rc.o
remove_antenna_LDFLAGS=remove_antenna.rc.o
fit_LDFLAGS=fit.rc.o
plot_ts_LDFLAGS=plot_ts.rc.o
SUFFIXES = .rc .o
RCFILES=bern2time.rc comfilt.rc comfilt_trend.rc fit.rc getpl.rc plot_ts.rc raw.rc remove_antenna.rc remove_period.rc remove_trend.rc

$(RCFILES): %.rc:
	printf "id ICON \"$$(basename $@ .rc).ico\"" > $@
	cat ../common/version.rc >> $@

bern2time.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

getpl.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

comfilt.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

raw.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

comfilt_trend.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

remove_trend.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

remove_period.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

remove_antenna.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

fit.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

plot_ts.rc.o: %.rc.o: %.rc
	$(WINDRES) $< -O coff -o $@

ico: %.rc
	convert -density 384 ../icons/sketchy_plot.svg -define icon:auto-resize $$(basename $< .rc).ico

clean-local:
	rm -f *.rc

endif
