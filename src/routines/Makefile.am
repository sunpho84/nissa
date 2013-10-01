AM_CPPFLAGS=-I$(top_srcdir)/src

#libraries
noinst_LTLIBRARIES=libroutines.la

#sources
libroutines_la_SOURCES= \
        ios.cpp \
        math_routines.cpp \
        mpi_routines.cpp

#headers
routinesincludedir=$(includedir)/routines
routinesinclude_HEADERS= \
	ios.hpp \
        math_routines.hpp \
        mpi_routines.hpp

#compile or not thread support
if USE_THREADS
libroutines_la_SOURCES+=thread.cpp
routinesinclude_HEADERS+=thread.hpp
endif
