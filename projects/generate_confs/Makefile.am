AM_CPPFLAGS=-I$(top_srcdir)/src

AUTOMAKE_OPTIONS=subdir-objects

LDADD=$(top_builddir)/lib/libnissa.a

bin_PROGRAMS= \
	$(top_builddir)/bin/generate_stag \
	$(top_builddir)/bin/generate_quenched


__top_builddir__bin_generate_stag_SOURCES=generate_stag.cpp
__top_builddir__bin_generate_stag_LDADD=$(LDADD)
__top_builddir__bin_generate_stag_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_generate_quenched_SOURCES=generate_quenched.cpp
__top_builddir__bin_generate_quenched_LDADD=$(LDADD)
__top_builddir__bin_generate_quenched_LDFLAGS=$(LDFLAGS)
