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
	tokenize_corr_fixed.cpp \
	parse_corr_fixed.cpp

parse_corr_fixed.cpp: $(top_srcdir)/config/fix_bison.sh Makefile parse_corr.cpp
	@echo "Fixing the parser"
	@bash $(top_srcdir)/config/fix_bison.sh parse_corr.cpp > $@
	@cp parse_corr_fixed.cpp parse_corr.cpp
	@touch parser_fixed.cpp parser.hpp

tokenize_corr_fixed.cpp: $(top_srcdir)/config/fix_flex.sh Makefile tokenize_corr.cpp
	@echo "Fixing the parser"
	@bash $(top_srcdir)/config/fix_flex.sh tokenize_corr.cpp > $@
	@cp tokenize_corr_fixed.cpp tokenize_corr.cpp
	@touch tokenizer_fixed.cpp

AM_YFLAGS=-d
BUILT_SOURCES=$(semileptonic_smeared_parser_built_sources)
CLEANFILES=$(semileptonic_smeared_parser_built_sources) 

###########################################################################

 __top_builddir__bin_semileptonic_smeared_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_sources)
 __top_builddir__bin_semileptonic_smeared_point_source_SOURCES=semileptonic_smeared.cpp $(semileptonic_smeared_parser_sources)
__top_builddir__bin_semileptonic_smeared_point_source_CXXFLAGS=-DPOINT_SOURCE_VERSION
