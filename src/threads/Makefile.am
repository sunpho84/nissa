################################### threads ###################################

include_HEADERS+= \
	%D%/cuda_threads.hpp \
	%D%/no_threads.hpp \
	%D%/openmp_threads.hpp \
	%D%/threads.hpp

if USE_OPENMP_THREADS
__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/openmp_threads.cpp
endif
