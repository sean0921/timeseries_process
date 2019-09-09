## GNU makefile

FC              = @FC@
FCFLAGS         = @FCFLAGS@
OBJS		= remove_antenna

all: $(OBJS)

$(OBJS): %: %.for
	$(FC) $(FCFLAGS) -o $@ $^

clean-local:
	rm -f $(OBJS) *.exe
