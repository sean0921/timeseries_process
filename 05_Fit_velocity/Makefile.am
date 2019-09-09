## GNU makefile

FC              = @FC@
FCFLAGS         = @FCFLAGS@
OBJS		= fit

all: $(OBJS)

$(OBJS): %: %.for
	$(FC) $(FCFLAGS) -o $@ $^

clean-local:
	rm -f $(OBJS) *.exe
