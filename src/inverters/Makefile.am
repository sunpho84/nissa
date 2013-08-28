#compile or not bgq test and source
if BGQ
bgq_src_file= \
        twisted_mass/cgm_invert_tmQ2_bgq.cpp \
        twisted_mass/cgm_invert_tmQ2_bgq.h
else
bgq_src_file=
endif

noinst_LTLIBRARIES=libinverters.la
libinverters_la_SOURCES= \
        twisted_mass/cg_invert_tmDeoimpr.h \
        twisted_mass/cg_invert_tmQ2.h \
        twisted_mass/cg_64_invert_tmQ2.h \
        twisted_mass/cgm_invert_tmQ2.h \
        twisted_mass/cg_128_invert_tmQ2.h \
        twisted_mass/cg_128_invert_tmDeoimpr.h \
        twisted_mass/tm_frontends.h \
        staggered/cgm_invert_stD2ee_m2.h \
        staggered/cg_invert_stD2ee_m2.h \
        staggered/cg_invert_stD.h \
        staggered/cg_invert_evn_stD.h \
        Wclov/cg_invert_WclovQ.h \
        Wclov/cg_invert_WclovQ2.h \
        tmclov/cg_invert_tmclovQ.h \
        tmclov/cg_invert_tmclovQ2.h \
        tmclov/cgm_invert_tmclovQ2.h \
        tmclov/tmclov_frontends.h \
        twisted_mass/cg_invert_tmDeoimpr.cpp \
        twisted_mass/cg_invert_tmQ2.cpp \
        twisted_mass/cg_64_invert_tmQ2.cpp \
        twisted_mass/cgm_invert_tmQ2.cpp \
        twisted_mass/cg_128_invert_tmQ2.cpp \
        twisted_mass/cg_128_invert_tmDeoimpr.cpp \
        twisted_mass/tm_frontends.cpp \
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

EXTRA_DIST= \
        twisted_mass/cg_invert_tmQ2_common.cpp \
	twisted_mass/cgm_invert_tmQ2_bgq.h \
        templates/cg_128_invert_template.cpp \
        templates/cg_invert_template.cpp \
        templates/cgm_invert_template.cpp \
        templates/cg_128_invert_template_threaded.cpp \
        templates/cg_invert_template_threaded.cpp \
        templates/cgm_invert_template_threaded.cpp
