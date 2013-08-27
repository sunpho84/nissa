#compile or not thread support
if USE_THREADS
thread_src=thread.cpp thread.h
else
thread_src=
endif

noinst_LTLIBRARIES=libroutines.la
libroutines_la_SOURCES= \
        ios.h \
        math_routines.h \
        mpi_routines.h \
        ios.cpp \
        math_routines.cpp \
        mpi_routines.cpp \
        $(thread_src)
