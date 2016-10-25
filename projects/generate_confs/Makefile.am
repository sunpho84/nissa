AM_CPPFLAGS=-I$(top_srcdir)/src

AUTOMAKE_OPTIONS=subdir-objects

LDADD=$(top_builddir)/lib/libnissa.a

bin_PROGRAMS= \
	$(top_builddir)/bin/generate \
	$(top_builddir)/bin/generate_quenched \
	$(top_builddir)/bin/convert_input

############################################### parser ####################################

generate_parser_sources= \
        parser.ypp \
        tokenizer.lpp \
        driver.hpp \
        redefine_yy.hpp
generate_parser_built_sources= \
        tokenizer.cpp \
        parser.hpp \
        parser.cpp \
	parser.cpp \
	tokenizer.cpp

AM_YFLAGS=-d
tokenizer.$(OBJEXT) : CXXFLAGS+=-Wno-sign-compare
BUILT_SOURCES=$(generate_parser_built_sources)
CLEANFILES=$(generate_parser_built_sources)

#############################################################################################

__top_builddir__bin_generate_SOURCES=generate.cpp $(generate_parser_sources)
__top_builddir__bin_generate_LDADD=$(LDADD)
__top_builddir__bin_generate_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_generate_quenched_SOURCES=generate_quenched.cpp
__top_builddir__bin_generate_quenched_LDADD=$(LDADD)
__top_builddir__bin_generate_quenched_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_convert_input_SOURCES=convert_input.cpp $(generate_parser_sources)
__top_builddir__bin_convert_input_LDADD=$(LDADD)
__top_builddir__bin_convert_input_LDFLAGS=$(LDFLAGS)
