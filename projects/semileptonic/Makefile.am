AM_CPPFLAGS=-I$(top_srcdir)/src

LDADD=$(top_builddir)/lib/libnissa.a

AUTOMAKE_OPTIONS=subdir-objects

#target programs
bin_PROGRAMS= \
	$(top_builddir)/bin/semileptonic_smeared \
	$(top_builddir)/bin/semileptonic_smeared_point_source

#semileptonic_smeared
semileptonic_smeared_parser_sources= \
        parse_corr.ypp \
        tokenize_corr.lpp \
        parse_corr.hpp \
        driver_corr.hpp \
	redefine_yy.hpp
semileptonic_smeared_parser_built_sources= \
	tokenize_corr.cpp \
	parse_corr.hpp \
	parse_corr.cpp

#with or without parser
AM_YFLAGS=-d
if PARSER_GENERATION
 BUILT_SOURCES=$(semileptonic_smeared_parser_built_sources)
 __top_builddir__bin_semileptonic_smeared_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_sources)
 __top_builddir__bin_semileptonic_smeared_point_source_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_sources)
else
 __top_builddir__bin_semileptonic_smeared_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_built_sources)
 __top_builddir__bin_semileptonic_smeared_point_source_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_built_sources)
endif
__top_builddir__bin_semileptonic_smeared_point_source_CXXFLAGS=-DPOINT_SOURCE_VERSION
