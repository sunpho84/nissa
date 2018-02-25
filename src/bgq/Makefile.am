######################################### bgq specific #####################################
if BGQ
__top_builddir__lib_libnissa_a_SOURCES+= \
	bgq/Wilson_hopping_matrix_lx_bgq.cpp \
	bgq/Wilson_hopping_matrix_eo_or_oe_bgq.cpp \
	bgq/staggered_hopping_matrix_eo_or_oe_bgq.cpp

include_HEADERS+= \
	bgq/Wilson_hopping_matrix_lx_bgq.hpp \
	bgq/Wilson_hopping_matrix_eo_or_oe_bgq_32_64.hpp \
	bgq/Wilson_hopping_matrix_eo_or_oe_bgq.hpp \
	bgq/staggered_hopping_matrix_eo_or_oe_bgq_32_64.hpp \
	bgq/staggered_hopping_matrix_eo_or_oe_bgq.hpp \
	bgq/bgq_macros.hpp \
	bgq/intrinsic/declare.hpp \
	bgq/intrinsic/load.hpp \
	bgq/intrinsic/mergesplit.hpp \
	bgq/intrinsic/oper.hpp \
	bgq/intrinsic/store.hpp \
	bgq/intrinsic.hpp

#bgq barrier
if BGQ_EMU
else
__top_builddir__lib_libnissa_a_SOURCES+=bgq/bgq_barrier.c
include_HEADERS+=bgq/bgq_barrier.hpp
endif

#compile or not spi stuff
if SPI
__top_builddir__lib_libnissa_a_SOURCES+=bgq/spi.cpp
include_HEADERS+=bgq/spi.hpp
endif

endif
