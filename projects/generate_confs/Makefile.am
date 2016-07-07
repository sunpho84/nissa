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
        parser.cpp

parser_fixed.cpp: $(top_srcdir)/projects/generate_confs/fix_bison.sh Makefile $(top_builddir)/projects/generate_confs/parser.cpp
	@echo "Fixing the parser"
	@bash $(top_srcdir)/projects/generate_confs/fix_bison.sh parser.cpp > $@

tokenizer_fixed.cpp: $(top_srcdir)/projects/generate_confs/fix_flex.sh Makefile $(top_builddir)/projects/generate_confs/tokenizer.cpp
	@echo "Fixing the parser"
	@bash $(top_srcdir)/projects/generate_confs/fix_flex.sh tokenizer.cpp > $@

AM_YFLAGS=-d
BUILT_SOURCES=$(generate_parser_built_sources) parser_fixed.cpp tokenizer_fixed.cpp

CLEANFILES=$(generate_parser_built_sources) parser_fixed.cpp tokenizer_fized.cpp

#############################################################################################

__top_builddir__bin_generate_SOURCES=generate.cpp parser.hpp parser_fixed.cpp tokenizer_fixed.cpp
__top_builddir__bin_generate_LDADD=$(LDADD)
__top_builddir__bin_generate_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_generate_quenched_SOURCES=generate_quenched.cpp
__top_builddir__bin_generate_quenched_LDADD=$(LDADD)
__top_builddir__bin_generate_quenched_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_convert_input_SOURCES=convert_input.cpp $(generate_parser_sources)
__top_builddir__bin_convert_input_LDADD=$(LDADD)
__top_builddir__bin_convert_input_LDFLAGS=$(LDFLAGS)
