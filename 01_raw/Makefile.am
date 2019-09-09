## GNU makefile

FC		= @FC@
FCFLAGS 	= @FCFLAGS@

OBJS		= comfilt comfilt_trend

all: $(OBJS)

$(OBJS): %: %.for
	$(FC) $(FCFLAGS) -o $@ $^

clean-local:
	rm -f $(OBJS) *.exe
