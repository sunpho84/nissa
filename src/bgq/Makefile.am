AM_CXXFLAGS=-I$(top_srcdir)/src

#bgq barrier
if BGQ_EMU
else
bgq_barrier_src_files=bgq_barrier.c
bgq_barrier_header_files=bgq_barrier.h
endif

#compile or not spi stuff
if SPI
spi_src_files=spi.cpp
spi_header_files=spi.h
endif

noinst_LTLIBRARIES=libbgq.la
libbgq_la_SOURCES= \
	$(spi_src_files) \
	$(bgq_barrier_src_files) \
        Wilson_hopping_matrix_lx_bgq.cpp \
        staggered_hopping_matrix_eo_or_oe_bgq.cpp

bgqincludedir=$(includedir)/bgq
bgqinclude_HEADERS= \
	Wilson_hopping_matrix_lx_bgq.h \
        staggered_hopping_matrix_eo_or_oe_bgq.h \
        bgq_macros.h \
        intrinsic/declare.h \
        intrinsic/load.h \
        intrinsic/mergesplit.h \
        intrinsic/oper.h \
        intrinsic/store.h \
        intrinsic.h
