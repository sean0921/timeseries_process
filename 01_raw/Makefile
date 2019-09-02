## GNU makefile

FC		:= gfortran
FC_MINGW	:= i686-w64-mingw32.static-gfortran
FCFLAGS 	:= -std=legacy -g -static -cpp
FCFLAGS_MINGW 	:= -DMINGW

OBJS		= comfilt  plot_ts  raw  comfilt_trend
OBJS_MINGW	= comfilt.exe  plot_ts.exe  raw.exe  comfilt_trend.exe

PROJECT_PATH	= 01_raw

all: $(OBJS)

%: %.for
	$(FC) $(FCFLAGS) -o $@ $^

raw: comfilt.for
	$(FC) $(FCFLAGS) -DRAW -o $@ $^

plot_ts: ../common/plot_ts.for
	$(FC) $(FCFLAGS) -o $@ $^

mingw: $(OBJS_MINGW)

%.exe: %.for
	$(FC_MINGW) $(FCFLAGS) $(FCFLAGS_MINGW) -o $@ $^

raw.exe: comfilt.for
	$(FC_MINGW) $(FCFLAGS) $(FCFLAGS_MINGW) -DRAW -o $@ $^

plot_ts.exe: ../common/plot_ts.for
	$(FC_MINGW) $(FCFLAGS) $(FCFLAGS_MINGW) -o $@ $^

package: $(OBJS)
	install -d ../build/$(PROJECT_PATH)/
	install -m 0755 $(OBJS) ../build/$(PROJECT_PATH)/
	install -m 0644 *.example ../build/$(PROJECT_PATH)/

mingw_package: $(OBJS_MINGW)
	install -d ../build/$(PROJECT_PATH)/
	install -m 0755 $(OBJS_MINGW) ../build/$(PROJECT_PATH)/
	install -m 0644 *.example ../build/$(PROJECT_PATH)/

help:
	@echo "all:           compile program for Linux or other *nix-like platform"
	@echo "mingw:         compile program for Win32 platform"
	@echo "package:       package program for *nix-like platform to build path"
	@echo "mingw_package: package program for Win32 platform to build path"

clean:
	rm -f $(OBJS) *.exe
