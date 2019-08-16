FC		:= gfortran
FC_MINGW	:= i686-w64-mingw32.static-gfortran
FCFLAGS 	:= -std=legacy -g -static

OBJS		= comfilt  plot_ts
OBJS_MINGW	= comfilt.exe  plot_ts.exe

all: $(OBJS)

%: %.for
	$(FC) $(FCFLAGS) -o $@ $^

mingw: $(OBJS_MINGW)

%.exe: %
	$(FC_MINGW) $(FCFLAGS) -o $@ $^

clean:
	rm -f $(OBJS) *.exe

install:
	install -m 0755 $(OBJS) test/

mingw_install:
	install -m 0755 $(OBJS_MINGW) test/
