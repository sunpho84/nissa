AM_CPPFLAGS=-I$(top_srcdir)/src

AUTOMAKE_OPTIONS=subdir-objects

LDADD=$(top_builddir)/lib/libnissa.a

bin_PROGRAMS= \
	$(top_builddir)/bin/checks \
	$(top_builddir)/bin/disco \
	$(top_builddir)/bin/ib \
	$(top_builddir)/bin/test_prop

common_sources=conf.cpp contr.cpp pars.cpp prop.cpp

include_HEADERS=conf.hpp contr.hpp pars.hpp prop.hpp


__top_builddir__bin_ib_SOURCES=ib.cpp $(common_sources)
__top_builddir__bin_ib_LDADD=$(LDADD)
__top_builddir__bin_ib_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_disco_SOURCES=disco.cpp $(common_sources)
__top_builddir__bin_disco_LDADD=$(LDADD)
__top_builddir__bin_disco_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_test_prop_SOURCES=test_prop.cpp $(common_sources)
__top_builddir__bin_test_prop_LDADD=$(LDADD)
__top_builddir__bin_test_prop_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_checks_SOURCES=checks.cpp $(common_sources)
__top_builddir__bin_checks_LDADD=$(LDADD)
__top_builddir__bin_checks_LDFLAGS=$(LDFLAGS)
