AM_CPPFLAGS=-I$(top_srcdir)/src

AUTOMAKE_OPTIONS=subdir-objects

LDADD=$(top_builddir)/lib/libnissa.a

#target programs
bin_PROGRAMS= \
	$(top_builddir)/bin/RIMOM

__top_builddir__bin_RIMOM_SOURCES=RIMOM.cpp
__top_builddir__bin_RIMOM_LDADD=$(LDADD)
__top_builddir__bin_RIMOM_LDFLAGS=$(LDFLAGS)
