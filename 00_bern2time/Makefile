FC		:= gfortran
FC_MINGW	:= i686-w64-mingw32.static-gfortran

OBJS		= bern2time_FN  bern2time
OBJS_MINGW	= bern2time_FN.exe  bern2time.exe

all: $(OBJS)

%: %.for
	$(FC) -o $@ $^

mingw: $(OBJS_MINGW)

%.exe: %.for
	$(FC_MINGW) -o $@ $^

clean:
	rm -f $(OBJS) *.exe
