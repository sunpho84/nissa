include $(top_srcdir)/config/baseMakefile.am

%D%/gitInfo.hpp: $(top_srcdir)/.git/HEAD $(top_srcdir)/.git/index $(top_srcdir)/config/getGitInfo.sh Makefile
	@echo "Reading git version"
	@$(top_srcdir)/config/getGitInfo.sh $(top_srcdir) > $@
$(top_srcdir)/.git/HEAD:
$(top_srcdir)/.git/index:

BUILT_SOURCES+= \
	%D%/gitInfo.hpp

include_HEADERS+= \
	nissa.hpp \
	%D%/gitInfo.hpp

include $(top_srcdir)/src/base/Makefile.am
include $(top_srcdir)/src/communicate/Makefile.am
include $(top_srcdir)/src/io/Makefile.am
include $(top_srcdir)/src/metaprogramming/Makefile.am
include $(top_srcdir)/src/newTypes/Makefile.am
include $(top_srcdir)/src/operations/Makefile.am
include $(top_srcdir)/src/routines/Makefile.am
include $(top_srcdir)/src/threads/Makefile.am
include $(top_srcdir)/src/tuples/Makefile.am
