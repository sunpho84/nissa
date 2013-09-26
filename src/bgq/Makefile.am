AM_CPPFLAGS=-I$(top_srcdir)/src

#bgq barrier
if BGQ_EMU
else
bgq_barrier_src_files=bgq_barrier.c
bgq_barrier_header_files=bgq_barrier.hpp
endif

#compile or not spi stuff
if SPI
spi_src_files=spi.cpp
spi_header_files=spi.hpp
endif

noinst_LTLIBRARIES=libbgq.la
libbgq_la_SOURCES= \
	$(spi_src_files) \
	$(bgq_barrier_src_files) \
        Wilson_hopping_matrix_lx_bgq.cpp \
        staggered_hopping_matrix_eo_or_oe_bgq.cpp

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
