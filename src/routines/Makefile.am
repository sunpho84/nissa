AM_CXXFLAGS=-I$(top_srcdir)/src

#compile or not thread support
if USE_THREADS
thread_src_files=thread.cpp
thread_headers_files=thread.h
endif

noinst_LTLIBRARIES=libroutines.la
libroutines_la_SOURCES= \
        ios.cpp \
        math_routines.cpp \
        mpi_routines.cpp \
        $(thread_src_files)

routinesincludedir=$(includedir)/routines
routinesinclude_HEADERS= \
	ios.h \
        math_routines.h \
        mpi_routines.h \
	$(thread_headers_files)