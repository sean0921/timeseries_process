## GNU makefile

FC              = @FC@
FCFLAGS         = @FCFLAGS@
OBJS		= plot_ts

all: $(OBJS)

$(OBJS): %: %.for
	$(FC) $(FCFLAGS) -o $@ $^

clean-local:
	rm -f $(OBJS) *.exe
