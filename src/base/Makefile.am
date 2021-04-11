########################################### base sources ##################################

%D%/git_info.hpp: $(top_srcdir)/.git/HEAD $(top_srcdir)/.git/index $(top_srcdir)/config/get_git_info.sh Makefile
	@echo "Reading git version"
	@$(top_srcdir)/config/get_git_info.sh $(top_srcdir) > $@
$(top_srcdir)/.git/HEAD:
$(top_srcdir)/.git/index:


BUILT_SOURCES=%D%/git_info.hpp

__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/bench.cpp \
	%D%/close.cpp \
	%D%/debug.cpp \
	%D%/init.cpp \
	%D%/vectors.cpp
include_HEADERS+= \
	%D%/bench.hpp \
	%D%/close.hpp \
	%D%/debug.hpp \
	%D%/git_info.hpp \
	%D%/init.hpp \
	%D%/vectors.hpp

if USE_TMLQCD
__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/tmLQCD_bridge.cpp
include_HEADERS+= \
	%D%/tmLQCD_bridge.hpp
endif

if USE_DDALPHAAMG
__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/DDalphaAMG_bridge.cpp
include_HEADERS+= \
	%D%/DDalphaAMG_bridge.hpp
endif

if USE_CUDA
__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/cuda.cpp
include_HEADERS+= \
	%D%/cuda.hpp
endif

if USE_QUDA
__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/quda_bridge.cpp
include_HEADERS+= \
	%D%/quda_bridge.hpp
endif
