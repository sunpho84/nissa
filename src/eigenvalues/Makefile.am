###################################### eigenvalues #################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/eigenvalues_autarchic.cpp \
	%D%/eigenvalues_overlap.cpp \
	%D%/eigenvalues_staggered.cpp

include_HEADERS+= \
	%D%/eigenvalues.hpp \
	%D%/eigenvalues_all.hpp \
	%D%/eigenvalues_autarchic.hpp \
	%D%/eigenvalues_overlap.hpp \
	%D%/eigenvalues_staggered.hpp

if USE_PARPACK
__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/eigenvalues_parpack.cpp

include_HEADERS+= \
	%D%/eigenvalues_parpack.hpp
endif
