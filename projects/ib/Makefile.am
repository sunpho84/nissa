AM_CPPFLAGS=-I$(top_srcdir)/src

AUTOMAKE_OPTIONS=subdir-objects

LDADD=$(top_builddir)/lib/libnissa.a

bin_PROGRAMS= \
	$(top_builddir)/bin/mel \
	$(top_builddir)/bin/mel_point_source \
	$(top_builddir)/bin/test_prop

common_sources=conf.cpp

include_HEADERS=conf.hpp pars.hpp


__top_builddir__bin_mel_SOURCES=mel.cpp $(common_sources)
__top_builddir__bin_mel_LDADD=$(LDADD)
__top_builddir__bin_mel_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_mel_point_source_SOURCES=mel.cpp $(common_sources)
__top_builddir__bin_mel_point_source_LDADD=$(LDADD)
__top_builddir__bin_mel_point_source_CXXFLAGS="-DPOINT_SOURCE_VERSION"
__top_builddir__bin_mel_point_source_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_test_prop_SOURCES=test_prop.cpp $(common_sources)
__top_builddir__bin_test_prop_LDADD=$(LDADD)
__top_builddir__bin_test_prop_LDFLAGS=$(LDFLAGS)
