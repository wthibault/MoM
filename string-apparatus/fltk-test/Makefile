include Makefile.FLTK

RM    = rm -f
SHELL = /bin/sh
.SILENT:

# Executables
ALL = hello$(EXEEXT) \
	hecko$(EXEEXT)

# default target -- build everything
default all: $(ALL)

# clean everything
clean:
	$(RM) $(ALL)
	$(RM) *.o
	$(RM) core
	$(RM) *~

