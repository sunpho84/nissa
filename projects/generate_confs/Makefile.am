bin_PROGRAMS+= \
	$(top_builddir)/bin/generate \
	$(top_builddir)/bin/convert_input

#	$(top_builddir)/bin/generate_quenched

############################################### parser ####################################

generate_parser_sources= \
        %D%/parser.ypp \
        %D%/tokenizer.lpp \
        %D%/driver.hpp \
        %D%/redefine_yy.hpp
generate_parser_built_sources= \
        %D%/tokenizer.cpp \
        %D%/parser.hpp \
        %D%/parser.cpp \
	%D%/tokenizer.cpp

%C%_tokenizer.$(OBJEXT) : CXXFLAGS+=-Wno-sign-compare
BUILT_SOURCES+=$(generate_parser_built_sources)
CLEANFILES+=$(generate_parser_built_sources)

#############################################################################################

__top_builddir__bin_generate_SOURCES=%D%/generate.cpp $(generate_parser_sources)
__top_builddir__bin_generate_LDADD=$(LDADD)
__top_builddir__bin_generate_LDFLAGS=$(LDFLAGS)

# __top_builddir__bin_generate_quenched_SOURCES=%D%/generate_quenched.cpp
# __top_builddir__bin_generate_quenched_LDADD=$(LDADD)
# __top_builddir__bin_generate_quenched_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_convert_input_SOURCES=%D%/convert_input.cpp $(generate_parser_sources)
__top_builddir__bin_convert_input_LDADD=$(LDADD)
__top_builddir__bin_convert_input_LDFLAGS=$(LDFLAGS)
