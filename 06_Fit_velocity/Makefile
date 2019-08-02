FC		:= gfortran
FC_MINGW	:= i686-w64-mingw32.static-gfortran

OBJS		= plot_ts  fit
OBJS_MINGW	= plot_ts.exe  fit.exe

all: $(OBJS)

%: %.for
	$(FC) -o $@ $^

mingw: $(OBJS_MINGW)

%.exe: %.for
	$(FC_MINGW) -o $@ $^

clean:
	rm -f $(OBJS) *.exe
