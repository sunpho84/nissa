AM_CPPFLAGS=-I$(top_srcdir)/src

#library naming
noinst_LTLIBRARIES=libbgq.la

#sources
libbgq_la_SOURCES= \
        Wilson_hopping_matrix_lx_bgq.cpp \
        staggered_hopping_matrix_eo_or_oe_bgq.cpp

#headers
bgqincludedir=$(includedir)/bgq
bgqinclude_HEADERS= \
	Wilson_hopping_matrix_lx_bgq.hpp \
        staggered_hopping_matrix_eo_or_oe_bgq.hpp \
        bgq_macros.hpp \
        intrinsic/declare.hpp \
        intrinsic/load.hpp \
        intrinsic/mergesplit.hpp \
        intrinsic/oper.hpp \
        intrinsic/store.hpp \
        intrinsic.hpp

#bgq barrier
if BGQ_EMU
else
libbgq_la_SOURCES+=bgq_barrier.c
bgqinclude_HEADERS+==bgq_barrier.hpp
endif

#compile or not spi stuff
if SPI
libbgq_la_SOURCES+=spi.cpp
bgqinclude_HEADERS+=spi.hpp
endif
