## GNU makefile

FC		= @FC@
FCFLAGS 	= @FCFLAGS@

OBJS		= bern2time  getpl

all: $(OBJS)

$(OBJS): %: %.for
	$(FC) $(FCFLAGS) -o $@ $^

clean-local:
	rm -f $(OBJS) *.exe
