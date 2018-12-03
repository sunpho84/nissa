###################################### eigenvalues #################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/eigenvalues_autarchic.cpp

include_HEADERS+= \
	%D%/eigenvalues.hpp \
	%D%/eigenvalues_autarchic.hpp

if USE_PARPACK
__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/eigenvalues_parpack.cpp

include_HEADERS+= \
	%D%/eigenvalues_parpack.hpp
endif
