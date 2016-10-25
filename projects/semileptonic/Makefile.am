AM_CPPFLAGS=-I$(top_srcdir)/src

AUTOMAKE_OPTIONS=subdir-objects

LDADD=$(top_builddir)/lib/libnissa.a

bin_PROGRAMS= \
	$(top_builddir)/bin/semileptonic_smeared \
	$(top_builddir)/bin/semileptonic_smeared_point_source

################################## parser ################################

semileptonic_smeared_parser_sources= \
        parse_corr.ypp \
        tokenize_corr.lpp \
        parse_corr.hpp \
        driver_corr.hpp \
	redefine_yy.hpp
semileptonic_smeared_parser_built_sources= \
	tokenize_corr.cpp \
	parse_corr.hpp \
	parse_corr.cpp \
	tokenize_corr.cpp \
	parse_corr.cpp

AM_YFLAGS=-d
tokenize_corr.$(OBJEXT) : CXXFLAGS+=-Wno-sign-compare
__top_builddir__bin_semileptonic_smeared_point_source-tokenize_corr.$(OBJEXT) : CXXFLAGS+=-Wno-sign-compare
BUILT_SOURCES=$(semileptonic_smeared_parser_built_sources)
CLEANFILES=$(semileptonic_smeared_parser_built_sources) 

###########################################################################

 __top_builddir__bin_semileptonic_smeared_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_sources)
 __top_builddir__bin_semileptonic_smeared_point_source_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_sources)
__top_builddir__bin_semileptonic_smeared_point_source_CXXFLAGS=-DPOINT_SOURCE_VERSION
