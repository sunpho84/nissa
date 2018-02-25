########################################## routines ###################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	routines/ios.cpp \
	routines/math_routines.cpp \
	routines/mpi_routines.cpp

include_HEADERS+= \
	routines/ios.hpp \
	routines/math_routines.hpp \
	routines/mpi_routines.hpp

#compile or not thread support
if USE_THREADS
__top_builddir__lib_libnissa_a_SOURCES+=routines/thread.cpp
include_HEADERS+=routines/thread.hpp
endif
