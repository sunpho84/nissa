bin_PROGRAMS+= \
	$(top_builddir)/bin/semileptonic_smeared \
	$(top_builddir)/bin/semileptonic_smeared_point_source

################################## parser ################################

semileptonic_smeared_parser_sources= \
	%D%/parse_corr.ypp \
        %D%/tokenize_corr.lpp \
        %D%/driver_corr.hpp \
	%D%/redefine_yy.hpp
semileptonic_smeared_parser_built_sources= \
	%D%/tokenize_corr.cpp \
	%D%/parse_corr.hpp \
	%D%/parse_corr.cpp \
	%D%/tokenize_corr.cpp

tokenize_corr.$(OBJEXT) : CXXFLAGS+=-Wno-sign-compare
__top_builddir__bin_semileptonic_smeared_point_source-tokenize_corr.$(OBJEXT) : CXXFLAGS+=-Wno-sign-compare
BUILT_SOURCES+=$(semileptonic_smeared_parser_built_sources)
CLEANFILES+=$(semileptonic_smeared_parser_built_sources) 

###########################################################################

 __top_builddir__bin_semileptonic_smeared_SOURCES=%D%/semileptonic_smeared.cpp $(semileptonic_smeared_parser_sources)
 __top_builddir__bin_semileptonic_smeared_point_source_SOURCES=%D%/semileptonic_smeared.cpp $(semileptonic_smeared_parser_sources)
__top_builddir__bin_semileptonic_smeared_point_source_CXXFLAGS=-DPOINT_SOURCE_VERSION
