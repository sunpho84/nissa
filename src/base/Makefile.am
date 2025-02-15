########################################### base sources ##################################

#add nissa library
lib_LIBRARIES=$(top_builddir)/lib/libnissa.a

%D%/git_info.hpp: $(top_srcdir)/.git/HEAD $(top_srcdir)/.git/index $(top_srcdir)/config/get_git_info.sh Makefile
	@echo "Reading git version"
	@$(top_srcdir)/config/get_git_info.sh $(top_srcdir) > $@
$(top_srcdir)/.git/HEAD:
$(top_srcdir)/.git/index:


BUILT_SOURCES+= \
	%D%/git_info.hpp

__top_builddir__lib_libnissa_a_SOURCES= \
	%D%/bench.cpp \
	%D%/debug.cpp \
	%D%/field.cpp \
	%D%/init.cpp \
	%D%/memory_manager.cpp \
	%D%/multiGridParams.cpp \
	%D%/random.cpp \
	%D%/vectors.cpp

include_HEADERS+= \
	%D%/bench.hpp \
	%D%/close.hpp \
	%D%/cuda_test.hpp \
	%D%/debug.hpp \
	%D%/field.hpp \
	%D%/git_info.hpp \
	%D%/init.hpp \
	%D%/memory_manager.hpp \
	%D%/multiGridParams.hpp \
	%D%/random.hpp \
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
include_HEADERS+= \
	%D%/cuda.hpp
endif

if USE_QUDA
__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/quda_bridge.cpp
include_HEADERS+= \
	%D%/quda_bridge.hpp
endif

if USE_EXTERNAL_SOLVER
__top_builddir__lib_libnissa_a_SOURCES+= \
	%D%/export_conf_to_external_solver.cpp
include_HEADERS+= \
	%D%/export_conf_to_external_solver.hpp
endif

EXTRA_DIST= \
	base/git_info.hpp
