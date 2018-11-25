######################################### bgq specific #####################################
if BGQ
__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/Wilson_hopping_matrix_lx_bgq.cpp \
	%D%/Wilson_hopping_matrix_eo_or_oe_bgq.cpp \
	%D%/staggered_hopping_matrix_eo_or_oe_bgq.cpp

include_HEADERS+= \
	%D%/Wilson_hopping_matrix_lx_bgq.hpp \
	%D%/Wilson_hopping_matrix_eo_or_oe_bgq_32_64.hpp \
	%D%/Wilson_hopping_matrix_eo_or_oe_bgq.hpp \
	%D%/staggered_hopping_matrix_eo_or_oe_bgq_32_64.hpp \
	%D%/staggered_hopping_matrix_eo_or_oe_bgq.hpp \
	%D%/bgq_macros.hpp \
	%D%/intrinsic/declare.hpp \
	%D%/intrinsic/load.hpp \
	%D%/intrinsic/mergesplit.hpp \
	%D%/intrinsic/oper.hpp \
	%D%/intrinsic/store.hpp \
	%D%/intrinsic.hpp

#bgq barrier
if BGQ_EMU
else
__top_builddir__lib_libnissa_a_SOURCES+=%D%/bgq_barrier.c
include_HEADERS+=%D%/bgq_barrier.hpp
endif

#compile or not spi stuff
if SPI
__top_builddir__lib_libnissa_a_SOURCES+=%D%/spi.cpp
include_HEADERS+=%D%/spi.hpp
endif

endif
