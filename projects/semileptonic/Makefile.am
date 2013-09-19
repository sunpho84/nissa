AM_CXXFLAGS=-I$(top_srcdir)/src

LDADD=$(top_builddir)/lib/libnissa.la

AUTOMAKE_OPTIONS=subdir-objects

#test bgq
if BGQ
test_bgq_target=test_bgq
endif

bin_PROGRAMS=$(top_builddir)/bin/semileptonic_smeared $(top_builddir)/bin/gevp $(test_bgq_target)

#include sources of test_bgq
if BGQ
test_bgq_SOURCES=test_bgq.cpp
endif

#semileptonic_smeared
semileptonic_smeared_parser_sources= \
        parse_corr.ypp \
        tokenize_corr.lpp \
        parse_corr.hpp \
        driver_corr.hpp \
	redefine_yy.h
semileptonic_smeared_parser_built_sources= \
	tokenize_corr.cpp \
	parse_corr.hpp \
	parse_corr.cpp
#with or without parser
AM_YFLAGS=-d
if PARSER_GENERATION
BUILT_SOURCES=$(semileptonic_smeared_parser_built_sources)
__top_builddir__bin_semileptonic_smeared_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_sources)
__top_builddir__bin_gevp_SOURCES=gevp.cpp $(semileptonic_smeared_parser_sources)
else
__top_builddir__bin_semileptonic_smeared_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_built_sources)
__top_builddir__bin_gevp_SOURCES=gevp.cpp $(semileptonic_smeared_parser_built_sources)
endif

#comment because it could happen that we don't have flex or bison
#CLEANFILES=$(semileptonic_smeared_parser_built_sources)
