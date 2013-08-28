LDADD=$(top_builddir)/lib/libnissa.la

AUTOMAKE_OPTIONS=subdir-objects
AM_YFLAGS=-d

semileptonic_smeared_parser_sources= \
        parse_corr.ypp \
        tokenize_corr.lpp \
        parse_corr.hpp \
        driver_corr.hpp
semileptonic_smeared_parser_built_sources= \
	tokenize_corr.cpp \
	parse_corr.hpp \
	parse_corr.cpp

BUILT_SOURCES=$(semileptonic_smeared_parser_built_sources)

bin_PROGRAMS=$(top_builddir)/bin/semileptonic_smeared

__top_builddir__bin_semileptonic_smeared_SOURCES= \
	semileptonic_smeared.cpp \
	$(semileptonic_smeared_parser_sources)

CLEANFILES=$(semileptonic_smeared_parser_built_sources)
