###################################### inverters #################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	inverters/momenta/cg_invert_MFACC.cpp \
	inverters/momenta/cgm_invert_MFACC.cpp \
	inverters/staggered/cgm_invert_stD2ee_m2_portable.cpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2.cpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2_portable.cpp \
	inverters/staggered/cg_invert_stD2ee_m2_portable.cpp \
	inverters/staggered/cg_invert_stD2Leb_ee_m2_portable.cpp \
	inverters/staggered/cgm_invert_stD2ee_m2.cpp \
	inverters/staggered/cg_invert_stD2ee_m2.cpp \
	inverters/staggered/cg_invert_stD.cpp \
	inverters/staggered/cg_invert_evn_stD.cpp \
	inverters/templates/eigenvalues.cpp \
	inverters/twisted_mass/cg_invert_tmD_eoprec.cpp \
	inverters/twisted_mass/cg_invert_tmQ2.cpp \
	inverters/twisted_mass/cg_64_invert_tmQ2.cpp \
	inverters/twisted_mass/cg_64_invert_tmD_eoprec.cpp \
	inverters/twisted_mass/cgm_invert_tmQ2.cpp \
	inverters/twisted_mass/cg_128_invert_tmQ2.cpp \
	inverters/twisted_mass/cg_128_invert_tmD_eoprec.cpp \
	inverters/twisted_mass/tm_frontends.cpp \
	inverters/twisted_clover/cg_invert_tmclovQ2.cpp \
	inverters/twisted_clover/cg_64_invert_tmclovQ2_portable.cpp \
	inverters/twisted_clover/cg_64_invert_tmclovD_eoprec.cpp \
	inverters/twisted_clover/cg_128_invert_tmclovQ2.cpp \
	inverters/twisted_clover/cg_128_invert_tmclovD_eoprec.cpp \
	inverters/twisted_clover/cgm_invert_tmclovDkern_eoprec_square_portable.cpp \
	inverters/twisted_clover/cgm_invert_tmclovQ2.cpp \
	inverters/twisted_clover/tmclov_frontends.cpp \
	inverters/Wclov/cg_invert_WclovQ.cpp \
	inverters/Wclov/cg_invert_WclovQ2.cpp \
	inverters/twisted_clover/cg_invert_tmclovD_eoprec.cpp \
	inverters/twisted_clover/cg_invert_tmclovQ.cpp

include_HEADERS+= \
	inverters/momenta/cg_invert_MFACC.hpp \
	inverters/momenta/cgm_invert_MFACC.hpp \
	inverters/templates/eigenvalues.hpp \
	inverters/twisted_mass/cg_invert_tmD_eoprec.hpp \
	inverters/twisted_mass/cg_invert_tmQ2.hpp \
	inverters/twisted_mass/cg_64_invert_tmQ2.hpp \
	inverters/twisted_mass/cg_64_invert_tmD_eoprec.hpp \
	inverters/twisted_mass/cgm_invert_tmQ2.hpp \
	inverters/twisted_mass/cg_128_invert_tmQ2.hpp \
	inverters/twisted_mass/cg_128_invert_tmD_eoprec.hpp \
	inverters/twisted_mass/tm_frontends.hpp \
	inverters/staggered/cgm_invert_stD2ee_m2.hpp \
	inverters/staggered/cg_invert_stD2ee_m2.hpp \
	inverters/staggered/cg_invert_stD.hpp \
	inverters/staggered/cg_invert_evn_stD.hpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2.hpp \
	inverters/staggered/cgm_invert_stD2ee_m2_portable.hpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2_portable.hpp \
	inverters/staggered/cg_invert_stD2ee_m2_portable.hpp \
	inverters/staggered/cg_invert_stD2Leb_ee_m2_portable.hpp \
	inverters/Wclov/cg_invert_WclovQ.hpp \
	inverters/Wclov/cg_invert_WclovQ2.hpp \
	inverters/twisted_clover/cg_invert_tmclovQ.hpp \
	inverters/twisted_clover/cg_invert_tmclovQ2.hpp \
	inverters/twisted_clover/cg_64_invert_tmclovQ2_portable.hpp \
	inverters/twisted_clover/cg_64_invert_tmclovQ2.hpp \
	inverters/twisted_clover/cg_64_invert_tmclovD_eoprec.hpp \
	inverters/twisted_clover/cg_128_invert_tmclovQ2.hpp \
	inverters/twisted_clover/cg_128_invert_tmclovD_eoprec.hpp \
	inverters/twisted_clover/cgm_invert_tmclovQ2.hpp \
	inverters/twisted_clover/cg_invert_tmclovD_eoprec.hpp \
	inverters/twisted_clover/tmclov_frontends.hpp

#compile or not bgq source
if BGQ
__top_builddir__lib_libnissa_a_SOURCES+= \
	inverters/staggered/cgm_invert_stD2ee_m2_bgq.cpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2_bgq.cpp \
	inverters/staggered/cg_invert_stD2ee_m2_bgq.cpp \
	inverters/twisted_mass/cgm_invert_tmQ2_bgq.cpp \
	inverters/twisted_mass/cg_64_invert_tmQ2_bgq.cpp \
	inverters/twisted_mass/cg_128_invert_tmQ2_bgq.cpp \
	inverters/twisted_mass/cg_64_invert_tmD_eoprec_bgq.cpp \
	inverters/twisted_mass/cg_invert_tmQ2_bgq.cpp \
	inverters/twisted_clover/cgm_invert_tmclovQ2_bgq.cpp \
	inverters/twisted_clover/cg_64_invert_tmclovQ2_bgq.cpp \
	inverters/twisted_clover/cg_64_invert_tmclovD_eoprec_bgq.cpp \
	inverters/twisted_clover/cg_128_invert_tmclovQ2_bgq.cpp

include_HEADERS+= \
	inverters/staggered/cgm_invert_stD2ee_m2_bgq.hpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2_bgq.hpp \
	inverters/staggered/cg_invert_stD2ee_m2_bgq.hpp \
	inverters/twisted_mass/cgm_invert_tmQ2_bgq.hpp \
	inverters/twisted_mass/cg_64_invert_tmQ2_bgq.hpp \
	inverters/twisted_mass/cg_128_invert_tmQ2_bgq.hpp \
	inverters/twisted_mass/cg_64_invert_tmD_eoprec_bgq.hpp \
	inverters/twisted_mass/cg_invert_tmQ2_bgq.hpp \
	inverters/twisted_clover/cg_64_invert_tmclovQ2_bgq.hpp \
	inverters/twisted_clover/cg_64_invert_tmclovD_eoprec_bgq.hpp \
	inverters/twisted_clover/cg_128_invert_tmclovQ2_bgq.hpp \
	inverters/twisted_clover/cgm_invert_tmclovQ2_bgq.hpp \
	inverters/twisted_clover/cgm_invert_tmclovDkern_eoprec_square_portable.hpp
endif

EXTRA_DIST+= \
	inverters/twisted_mass/cg_invert_tmQ2_common.cpp \
	inverters/twisted_mass/cgm_invert_tmQ2_bgq.hpp \
	inverters/templates/cgm_32_invert_template_threaded.cpp \
	inverters/templates/cg_128_invert_template_threaded.cpp \
	inverters/templates/cgm_invert_template_summsol_threaded.hpp \
	inverters/templates/cg_invert_template_threaded.cpp \
	inverters/templates/bicgstab_invert_template_threaded.cpp \
	inverters/templates/cgm_invert_template_threaded.cpp

