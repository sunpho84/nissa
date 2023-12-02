bin_PROGRAMS+= \
	$(top_builddir)/bin/test_wip

assembly_reports+=\
        %D%/test_wip.s

__top_builddir__bin_test_wip_SOURCES=%D%/test_wip.cpp
__top_builddir__bin_test_wip_LDADD=$(LDADD)
__top_builddir__bin_test_wip_LDFLAGS=$(LDFLAGS)

#generate assembly report if asked
if ASSEMBLY_REPORT
BUILT_SOURCES+=$(assembly_reports)
CLEANFILES+=$(assembly_reports)
endif
