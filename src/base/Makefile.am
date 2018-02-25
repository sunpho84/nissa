########################################### base sources ##################################

#add nissa library
lib_LIBRARIES=$(top_builddir)/lib/libnissa.a

base/git_info.hpp: $(top_srcdir)/.git/HEAD $(top_srcdir)/.git/index $(top_srcdir)/config/get_git_info.sh Makefile
	@echo "Reading git version"
	@$(top_srcdir)/config/get_git_info.sh $(top_srcdir) > $@
$(top_srcdir)/.git/HEAD:
$(top_srcdir)/.git/index:


BUILT_SOURCES=base/git_info.hpp

__top_builddir__lib_libnissa_a_SOURCES= \
	base/bench.cpp \
	base/close.cpp \
	base/debug.cpp \
	base/init.cpp \
	base/random.cpp \
	base/vectors.cpp
include_HEADERS+= \
	base/bench.hpp \
	base/close.hpp \
	base/debug.hpp \
	base/git_info.hpp \
	base/init.hpp \
	base/random.hpp \
	base/thread_macros.hpp \
	base/vectors.hpp

if USE_TMLQCD
__top_builddir__lib_libnissa_a_SOURCES+= \
	base/tmLQCD_bridge.cpp
include_HEADERS+= \
	base/tmLQCD_bridge.hpp
endif

if USE_DDALPHAAMG
__top_builddir__lib_libnissa_a_SOURCES+= \
	base/DDalphaAMG_bridge.cpp
include_HEADERS+= \
	base/DDalphaAMG_bridge.hpp
endif

#add BGQ if needed
if BGQ
include_HEADERS+=base/bgq_instructions.hpp
endif
