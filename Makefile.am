#### Simple Makefile for build a pack

SUBDIR = 00_bern2time 01_raw 02_remove_trend 03_remove_period 04_remove_antenna 05_Fit_velocity util

all:
	for i in $(SUBDIR);do cd $$i;$(MAKE);cd -;done

clean-local:
	for i in $(SUBDIR);do cd $$i;$(MAKE) clean;cd -;done
