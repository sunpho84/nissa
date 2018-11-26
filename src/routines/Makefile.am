########################################## routines ###################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/ios.cpp \
	%D%/math_routines.cpp \
	%D%/mpi_routines.cpp

include_HEADERS+= \
	%D%/ios.hpp \
	%D%/math_routines.hpp \
	%D%/mpi_routines.hpp

#compile or not thread support
if USE_THREADS
__top_builddir__lib_libnissa_a_SOURCES+=%D%/thread.cpp
include_HEADERS+=%D%/thread.hpp
endif
