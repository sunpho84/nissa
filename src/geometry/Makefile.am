################################### geometry ###################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/geometry_eo.cpp \
	%D%/geometry_Leb.cpp \
	%D%/geometry_lx.cpp \
	%D%/geometry_mix.cpp

include_HEADERS+= \
	%D%/geometry_eo.hpp \
	%D%/geometry_Leb.hpp \
	%D%/geometry_lx.hpp \
	%D%/geometry_mix.hpp

#compile or not virtual parallelization
if USE_VNODES
__top_builddir__lib_libnissa_a_SOURCES+=%D%/geometry_vir.cpp
include_HEADERS+=%D%/geometry_vir.hpp
endif
