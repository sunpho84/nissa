AM_CPPFLAGS=-I$(top_srcdir)/src

AUTOMAKE_OPTIONS=subdir-objects

LDADD=$(top_builddir)/lib/libnissa.a

bin_PROGRAMS= \
	$(top_builddir)/bin/checks \
	$(top_builddir)/bin/disco \
	$(top_builddir)/bin/mel \
	$(top_builddir)/bin/mel_point_source \
	$(top_builddir)/bin/nucleons \
	$(top_builddir)/bin/test_prop

common_sources=conf.cpp contr.cpp pars.cpp prop.cpp

include_HEADERS=conf.hpp contr.hpp pars.hpp prop.hpp


__top_builddir__bin_mel_SOURCES=mel.cpp $(common_sources)
__top_builddir__bin_mel_LDADD=$(LDADD)
__top_builddir__bin_mel_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_disco_SOURCES=disco.cpp $(common_sources)
__top_builddir__bin_disco_LDADD=$(LDADD)
__top_builddir__bin_disco_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_mel_point_source_SOURCES=mel.cpp $(common_sources)
__top_builddir__bin_mel_point_source_LDADD=$(LDADD)
__top_builddir__bin_mel_point_source_CXXFLAGS="-DPOINT_SOURCE_VERSION"
__top_builddir__bin_mel_point_source_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_nucleons_SOURCES=nucleons.cpp $(common_sources)
__top_builddir__bin_nucleons_LDADD=$(LDADD)
__top_builddir__bin_nucleons_CXXFLAGS="-DPOINT_SOURCE_VERSION"
__top_builddir__bin_nucleons_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_test_prop_SOURCES=test_prop.cpp $(common_sources)
__top_builddir__bin_test_prop_LDADD=$(LDADD)
__top_builddir__bin_test_prop_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_checks_SOURCES=checks.cpp $(common_sources)
__top_builddir__bin_checks_LDADD=$(LDADD)
__top_builddir__bin_checks_LDFLAGS=$(LDFLAGS)
