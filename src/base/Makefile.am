########################################### base sources ##################################

%D%/git_info.hpp: $(top_srcdir)/.git/HEAD $(top_srcdir)/.git/index $(top_srcdir)/config/get_git_info.sh Makefile
	@echo "Reading git version"
	@$(top_srcdir)/config/get_git_info.sh $(top_srcdir) > $@
$(top_srcdir)/.git/HEAD:
$(top_srcdir)/.git/index:


BUILT_SOURCES+= \
	%D%/git_info.hpp

include_HEADERS+= \
	%D%/close.hpp \
	%D%/cuda_test.hpp \
	%D%/debug.hpp \
	%D%/git_info.hpp \
	%D%/init.hpp \
	%D%/memory_manager.hpp \
	%D%/multiGridParams.hpp

