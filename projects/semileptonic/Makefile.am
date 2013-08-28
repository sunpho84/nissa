LDADD=$(top_builddir)/lib/libnissa.la

AUTOMAKE_OPTIONS=subdir-objects

bin_PROGRAMS=$(top_builddir)/bin/semileptonic_smeared

semileptonic_smeared_parser_sources= \
        parse_corr.ypp \
        tokenize_corr.lpp \
        parse_corr.hpp \
        driver_corr.hpp
semileptonic_smeared_parser_built_sources= \
	tokenize_corr.cpp \
	parse_corr.hpp \
	parse_corr.cpp

AM_YFLAGS=-d

if PARSER_GENERATION
BUILT_SOURCES=$(semileptonic_smeared_parser_built_sources)
__top_builddir__bin_semileptonic_smeared_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_sources)
else
__top_builddir__bin_semileptonic_smeared_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_built_sources)
endif


#comment because it could happen that we don't have flex or bison
#CLEANFILES=$(semileptonic_smeared_parser_built_sources)
