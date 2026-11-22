if USE_PARSE_PACT
bin_PROGRAMS+= \
	$(top_builddir)/bin/NICE

check_PROGRAMS= \
	$(top_builddir)/bin/nice_validate

TESTS=$(check_PROGRAMS)

nice_common_sources= \
	%D%/nice.hpp \
	%D%/checkpointIdentity.hpp \
	%D%/runOptions.hpp \
	%D%/runPolicy.hpp \
	%D%/rng.hpp \
	%D%/rngLayout.hpp \
	%D%/configurationSources.cpp \
	%D%/configurationSources.hpp \
	%D%/gaugeRecipes.hpp \
	%D%/tensorComps.hpp \
	%D%/tensorValues.hpp \
	%D%/sField.hpp \
	%D%/lineValue.hpp \
	%D%/taskExecution.hpp \
	%D%/taskExecution.cpp \
	%D%/programExecution.hpp \
	%D%/corrParts.hpp \
	%D%/correlations.hpp \
	%D%/graphCompiler.hpp \
	%D%/nASTy.cpp \
	%D%/nASTy.hpp \
	%D%/nASTyHelp.inc \
	%D%/scriptSource.hpp \
	%D%/scriptCompletion.hpp \
	%D%/scriptCompletion.cpp \
	%D%/niceSetup.cpp \
	%D%/niceSetup.hpp \
	%D%/correlationTaskGraph.hpp \
	%D%/compiledCorrelations.hpp \
	%D%/schedulingStrategies.hpp \
	%D%/legacyInputWriter.cpp \
	%D%/legacyInputWriter.hpp \
	%D%/corrParts.cpp

__top_builddir__bin_NICE_SOURCES= \
	%D%/NICE.cpp \
	%D%/checkpoint.hpp \
	%D%/wallTime.hpp \
	$(nice_common_sources)

__top_builddir__bin_nice_validate_SOURCES= \
	%D%/validate.cpp \
	$(nice_common_sources)
endif

EXTRA_DIST+= \
	%D%/nissa-script-mode.el \
	%D%/nissa-script-check.el \
	%D%/nissa-script-completion.el \
	%D%/nissa-completion.bash
