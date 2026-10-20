if USE_PARSE_PACT
bin_PROGRAMS+= \
	$(top_builddir)/bin/nice

__top_builddir__bin_nice_SOURCES= \
	%D%/nice.cpp \
	%D%/corrParts.cpp \
	%D%/corrParts.hpp \
	%D%/correlations.hpp \
	%D%/graphCompiler.hpp \
	%D%/nASTy.hpp \
	%D%/niceSetup.cpp \
	%D%/niceSetup.hpp \
	%D%/correlationTaskGraph.hpp \
	%D%/compiledCorrelations.hpp \
	%D%/schedulingStrategies.hpp \
	%D%/legacyInputWriter.cpp \
	%D%/legacyInputWriter.hpp

if USE_PREPARE_INPUT
__top_builddir__bin_nice_SOURCES+=%D%/prepareInput.cpp
endif
endif
