################################### geometry ###################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	geometry/geometry_eo.cpp \
	geometry/geometry_Leb.cpp \
	geometry/geometry_lx.cpp \
	geometry/geometry_mix.cpp

include_HEADERS+= \
	geometry/geometry_eo.hpp \
	geometry/geometry_Leb.hpp \
	geometry/geometry_lx.hpp \
	geometry/geometry_mix.hpp

#compile or not virtual parallelization
if USE_VNODES
__top_builddir__lib_libnissa_a_SOURCES+=geometry/geometry_vir.cpp
include_HEADERS+=geometry/geometry_vir.hpp
endif
