FC		:= gfortran
FC_MINGW	:= i686-w64-mingw32.static-gfortran

OBJS		= comfilt  plot_ts
OBJS_MINGW	= comfilt.exe  plot_ts.exe

all: $(OBJS)

%: %.for
	$(FC) -o $@ $^

mingw: $(OBJS_MINGW)

%.exe: %
	$(FC_MINGW) -o $@ $^

clean:
	rm -f $(OBJS) *.exe
