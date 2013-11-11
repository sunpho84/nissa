AM_CPPFLAGS=-I$(top_srcdir)/src

#compile or not bgq test and source
if BGQ
bgq_src_file= \
        staggered/cgm_invert_stD2ee_m2_bgq.cpp \
	staggered/cgm_invert_stD2ee_m2_bgq.hpp \
        staggered/cg_invert_stD2ee_m2_bgq.cpp \
	staggered/cg_invert_stD2ee_m2_bgq.hpp \
        twisted_mass/cgm_invert_tmQ2_bgq.cpp \
        twisted_mass/cgm_invert_tmQ2_bgq.hpp \
        twisted_mass/cg_128_invert_tmQ2_bgq.cpp \
        twisted_mass/cg_128_invert_tmQ2_bgq.hpp
endif

noinst_LTLIBRARIES=libinverters.la
libinverters_la_SOURCES= \
	$(bgq_src_file) \
        twisted_mass/cg_invert_tmDeoimpr.cpp \
        twisted_mass/cg_invert_tmQ2.cpp \
        twisted_mass/cg_64_invert_tmQ2.cpp \
        twisted_mass/cgm_invert_tmQ2.cpp \
        twisted_mass/cg_128_invert_tmQ2.cpp \
        twisted_mass/cg_128_invert_tmDeoimpr.cpp \
        twisted_mass/tm_frontends.cpp \
        staggered/cgm_invert_stD2ee_m2_portable.cpp \
        staggered/cg_invert_stD2ee_m2_portable.cpp \
        staggered/cgm_invert_stD2ee_m2.cpp \
        staggered/cg_invert_stD2ee_m2.cpp \
        staggered/cg_invert_stD.cpp \
        staggered/cg_invert_evn_stD.cpp \
        Wclov/cg_invert_WclovQ.cpp \
        Wclov/cg_invert_WclovQ2.cpp \
        tmclov/cg_invert_tmclovQ.cpp \
        tmclov/cg_invert_tmclovQ2.cpp \
        tmclov/cgm_invert_tmclovQ2.cpp \
        tmclov/tmclov_frontends.cpp

inverters_twisted_mass_includedir=$(includedir)/inverters/twisted_mass
inverters_twisted_mass_include_HEADERS= \
        twisted_mass/cg_invert_tmDeoimpr.hpp \
        twisted_mass/cg_invert_tmQ2.hpp \
        twisted_mass/cg_64_invert_tmQ2.hpp \
        twisted_mass/cgm_invert_tmQ2.hpp \
        twisted_mass/cg_128_invert_tmQ2.hpp \
        twisted_mass/cg_128_invert_tmDeoimpr.hpp \
        twisted_mass/tm_frontends.hpp

inverters_staggered_includedir=$(includedir)/inverters/staggered
inverters_staggered_include_HEADERS= \
        staggered/cgm_invert_stD2ee_m2.hpp \
        staggered/cg_invert_stD2ee_m2.hpp \
        staggered/cg_invert_stD.hpp \
        staggered/cg_invert_evn_stD.hpp \
	staggered/cgm_invert_stD2ee_m2_portable.hpp \
	staggered/cg_invert_stD2ee_m2_portable.hpp

inverters_Wclov_includedir=$(includedir)/inverters/Wclov
inverters_Wclov_include_HEADERS= \
        Wclov/cg_invert_WclovQ.hpp \
        Wclov/cg_invert_WclovQ2.hpp

inverters_tmclov_includedir=$(includedir)/inverters/tmclov
inverters_tmclov_include_HEADERS= \
        tmclov/cg_invert_tmclovQ.hpp \
        tmclov/cg_invert_tmclovQ2.hpp \
        tmclov/cgm_invert_tmclovQ2.hpp \
        tmclov/tmclov_frontends.hpp

EXTRA_DIST= \
        twisted_mass/cg_invert_tmQ2_common.cpp \
	twisted_mass/cgm_invert_tmQ2_bgq.hpp \
        templates/cg_128_invert_template_threaded.cpp \
        templates/cg_invert_template_threaded.cpp \
        templates/cgm_invert_template_threaded.cpp
