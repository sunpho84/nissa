if USE_PARSE_PACT
bin_PROGRAMS+= \
	$(top_builddir)/bin/nice

check_PROGRAMS= \
	$(top_builddir)/bin/nice_validate

TESTS=$(check_PROGRAMS)

nice_common_sources= \
	%D%/nice.hpp \
	%D%/runOptions.hpp \
	%D%/configurationSources.cpp \
	%D%/configurationSources.hpp \
	%D%/gaugeRecipes.hpp \
	%D%/tensorComps.hpp \
	%D%/tensorValues.hpp \
	%D%/sField.hpp \
	%D%/lineValue.hpp \
	%D%/taskExecution.hpp \
	%D%/programExecution.hpp \
	%D%/corrParts.hpp \
	%D%/correlations.hpp \
	%D%/graphCompiler.hpp \
	%D%/nASTy.cpp \
	%D%/nASTy.hpp \
	%D%/niceSetup.cpp \
	%D%/niceSetup.hpp \
	%D%/correlationTaskGraph.hpp \
	%D%/compiledCorrelations.hpp \
	%D%/schedulingStrategies.hpp \
	%D%/legacyInputWriter.cpp \
	%D%/legacyInputWriter.hpp \
	%D%/corrParts.cpp

__top_builddir__bin_nice_SOURCES= \
	%D%/nice.cpp \
	$(nice_common_sources)

__top_builddir__bin_nice_validate_SOURCES= \
	%D%/validate.cpp \
	$(nice_common_sources)
endif
