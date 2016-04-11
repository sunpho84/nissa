AM_CPPFLAGS=-I$(top_srcdir)/src

AUTOMAKE_OPTIONS=subdir-objects

LDADD=$(top_builddir)/lib/libnissa.a

bin_PROGRAMS= \
	$(top_builddir)/bin/nucleons_smeared

__top_builddir__bin_nucleons_smeared_SOURCES=nucleons_smeared.cpp
__top_builddir__bin_nucleons_smeared_LDADD=$(LDADD)
__top_builddir__bin_nucleons_smeared_LDFLAGS=$(LDFLAGS)
