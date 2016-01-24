AM_CPPFLAGS=-I$(top_srcdir)/src

AUTOMAKE_OPTIONS=subdir-objects

LDADD=$(top_builddir)/lib/libnissa.a

bin_PROGRAMS= \
	$(top_builddir)/bin/test_parser

#	$(top_builddir)/bin/generate_quenched
#	$(top_builddir)/bin/generate

############################################### parser ####################################

generate_parser_sources= \
        parser.ypp \
        tokenizer.lpp \
        parser.hpp \
        driver.hpp \
        redefine_yy.hpp
generate_parser_built_sources= \
        tokenizer.cpp \
        parser.hpp \
        parser.cpp

AM_YFLAGS=-d
BUILT_SOURCES=$(generate_parser_built_sources)

CLEANFILES=$(generate_parser_built_sources)

#############################################################################################

#__top_builddir__bin_generate_SOURCES=generate_.cpp $(generate_parser_sources)
#__top_builddir__bin_generate_LDADD=$(LDADD)
#__top_builddir__bin_generate_LDFLAGS=$(LDFLAGS)

#__top_builddir__bin_generate_quenched_SOURCES=generate_quenched.cpp
#__top_builddir__bin_generate_quenched_LDADD=$(LDADD)
#__top_builddir__bin_generate_quenched_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_test_parser_SOURCES=test_parser.cpp $(generate_parser_sources)
__top_builddir__bin_test_parser_LDADD=$(LDADD)
__top_builddir__bin_test_parser_LDFLAGS=$(LDFLAGS)
