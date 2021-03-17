bin_PROGRAMS+= \
	$(top_builddir)/bin/ib \
	$(top_builddir)/bin/disco

# $(top_builddir)/bin/checks_ib
# $(top_builddir)/bin/EU3_source
# 	$(top_builddir)/bin/disco3
#       $(top_builddir)/bin/mel_point_source
#	$(top_builddir)/bin/mel
#	$(top_builddir)/bin/disco
#	$(top_builddir)/bin/test_prop
#	$(top_builddir)/bin/checks
#	$(top_builddir)/bin/check_fftw

common_sources= \
	%D%/conf.cpp \
	%D%/contr.cpp \
	%D%/meslep.cpp \
	%D%/pars.cpp \
	%D%/prop.cpp

include_HEADERS+= \
	%D%/conf.hpp \
	%D%/contr.hpp \
	%D%/meslep.hpp \
	%D%/pars.hpp \
	%D%/prop.hpp

if USE_QUDA
bin_PROGRAMS+= \
	$(top_builddir)/bin/test_quda
 __top_builddir__bin_test_quda_SOURCES=%D%/test_quda.cpp
 __top_builddir__bin_test_quda_LDADD=$(LDADD)
 __top_builddir__bin_test_quda_LDFLAGS=$(LDFLAGS)
endif

 # __top_builddir__bin_checks_ib_SOURCES=%D%/checks.cpp
 # __top_builddir__bin_checks_ib_LDADD=$(LDADD)
 # __top_builddir__bin_checks_ib_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_ib_SOURCES=%D%/ib.cpp $(common_sources)
__top_builddir__bin_ib_LDADD=$(LDADD)
__top_builddir__bin_ib_LDFLAGS=$(LDFLAGS)

# __top_builddir__bin_EU3_source_SOURCES=%D%/EU3_source.cpp $(common_sources)
# __top_builddir__bin_EU3_source_LDADD=$(LDADD)
# __top_builddir__bin_EU3_source_LDFLAGS=$(LDFLAGS)

# __top_builddir__bin_mel_SOURCES=%D%/mel.cpp
# __top_builddir__bin_mel_LDADD=$(LDADD)
# __top_builddir__bin_mel_LDFLAGS=$(LDFLAGS)

# __top_builddir__bin_mel_point_source_SOURCES=%D%/mel.cpp
# __top_builddir__bin_mel_point_source_LDADD=$(LDADD)
# __top_builddir__bin_mel_point_source_CXXFLAGS="-DPOINT_SOURCE_VERSION"
# __top_builddir__bin_mel_point_source_LDFLAGS=$(LDFLAGS)

__top_builddir__bin_disco_SOURCES=%D%/disco.cpp $(common_sources)
__top_builddir__bin_disco_LDADD=$(LDADD)
__top_builddir__bin_disco_LDFLAGS=$(LDFLAGS)

# __top_builddir__bin_disco2_SOURCES=%D%/disco2.cpp
# __top_builddir__bin_disco2_LDADD=$(LDADD)
# __top_builddir__bin_disco2_LDFLAGS=$(LDFLAGS)

# __top_builddir__bin_disco3_SOURCES=%D%/disco3.cpp
# __top_builddir__bin_disco3_LDADD=$(LDADD)
# __top_builddir__bin_disco3_LDFLAGS=$(LDFLAGS)

# __top_builddir__bin_test_prop_SOURCES=%D%/test_prop.cpp $(common_sources)
# __top_builddir__bin_test_prop_LDADD=$(LDADD)
# __top_builddir__bin_test_prop_LDFLAGS=$(LDFLAGS)

# __top_builddir__bin_check_fftw_SOURCES=%D%/check_fftw.cpp
# __top_builddir__bin_check_fftw_ib_LDADD=$(LDADD)
# __top_builddir__bin_check_fftw_LDFLAGS=$(LDFLAGS)
