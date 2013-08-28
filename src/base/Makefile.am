BUILT_SOURCES=update_svnversion svnversion.h

if BGQ
bgq_src_files= \
	../bgq/spi.h
else
bgq_src_files=
endif


noinst_LTLIBRARIES=libbase.la
libbase_la_SOURCES= \
	bgq_instructions.h \
	close.cpp \
	close.h \
	debug.cpp \
	debug.h \
	global_variables.cpp \
	global_variables.h \
	init.cpp \
	init.h \
	macros.h \
	random.cpp \
	random.h \
	thread_macros.h \
	vectors.cpp \
	svnversion.h \
	vectors.h

EXTRA_DIST=update_svnversion.h.sh svnversion.h

update_svnversion:
	./update_svnversion.h.sh

.PHONY: update_svnversion