AM_CPPFLAGS=-I$(top_srcdir)/src

AUTOMAKE_OPTIONS=subdir-objects

LDADD=$(top_builddir)/lib/libnissa.a

bin_PROGRAMS= \
	$(top_builddir)/bin/mel \
	$(top_builddir)/bin/mel_point_source

__top_builddir__bin_mel_SOURCES=mel.cpp
__top_builddir__bin_mel_LDADD=$(LDADD)
__top_builddir__bin_mel_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_mel_point_source_SOURCES=mel.cpp
__top_builddir__bin_mel_point_source_LDADD=$(LDADD)
__top_builddir__bin_mel_point_source_CXXFLAGS="-DPOINT_SOURCE_VERSION"
__top_builddir__bin_mel_point_source_LDFLAGS=$(LDFLAGS)
