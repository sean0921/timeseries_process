## GNU makefile

FC		:= gfortran
FC_MINGW	:= i686-w64-mingw32.static-gfortran
FCFLAGS 	:= -std=legacy -g -static -cpp
FCFLAGS_MINGW 	:= -DMINGW

OBJS		= remove_trend  plot_ts
OBJS_MINGW	= remove_trend.exe  plot_ts.exe

all: $(OBJS)

%: %.for
	$(FC) $(FCFLAGS) -o $@ $^

mingw: $(OBJS_MINGW)

%.exe: %.for
	$(FC_MINGW) $(FCFLAGS) $(FCFLAGS_MINGW) -o $@ $^

install: $(OBJS)
	install -m 0755 $(OBJS) test/

mingw_install: $(OBJS_MINGW)
	install -m 0755 $(OBJS_MINGW) test/

help:
	@echo "all:           compile program for Linux or other *nix-like platform"
	@echo "mingw:         compile program for Win32 platform"
	@echo "install:       install program for *nix-like platform to test/"
	@echo "mingw_install: install program for Win32 platform to test/"

clean:
	rm -f $(OBJS) *.exe
