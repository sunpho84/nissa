######################################## dirac operators#####################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	dirac_operators/momenta/MFACC.cpp \
	dirac_operators/stD/dirac_operator_stD.cpp \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.cpp \
	dirac_operators/tmQ/dirac_operator_tmQ.cpp \
	dirac_operators/tmQ_left/dirac_operator_tmQ_left.cpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2.cpp \
	dirac_operators/tmQ/reconstruct_tm_doublet.cpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.cpp \
	dirac_operators/tmclovQ/reconstruct_tmclov_doublet.cpp \
	dirac_operators/Wstat/dirac_operator_Wstat.cpp \
	dirac_operators/WclovQ/dirac_operator_WclovQ.cpp \
	dirac_operators/WclovQ2/dirac_operator_WclovQ2.cpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ.cpp \
        dirac_operators/overlap/dirac_H_overlap.cpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ_128.cpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.cpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_128.cpp \
	dirac_operators/tmQ/dirac_operator_tmQ_128.cpp \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_128.cpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_128.cpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_128.cpp

include_HEADERS+= \
	dirac_operators/momenta/MFACC.hpp \
	dirac_operators/stD/dirac_operator_stD.hpp\
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_128.hpp \
	dirac_operators/tmQ/dirac_operator_tmQ.hpp \
	dirac_operators/tmQ/reconstruct_tm_doublet.hpp \
	dirac_operators/tmQ/dirac_operator_tmQ_128.hpp \
	dirac_operators/tmQ_left/dirac_operator_tmQ_left.hpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2.hpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_128.hpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ_128.hpp \
	dirac_operators/tmclovQ/reconstruct_tmclov_doublet.hpp \
	dirac_operators/Wstat/dirac_operator_Wstat.hpp \
	dirac_operators/WclovQ/dirac_operator_WclovQ.hpp \
	dirac_operators/WclovQ2/dirac_operator_WclovQ2.hpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.hpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_128.hpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_128.hpp \
        dirac_operators/overlap/dirac_H_overlap.hpp 
#compile or not bgq source
if BGQ
__top_builddir__lib_libnissa_a_SOURCES+= \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_bgq.cpp \
	dirac_operators/tmQ/dirac_operator_tmQ_bgq.cpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_bgq.cpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_128_bgq.cpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_bgq.cpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ_bgq.cpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_bgq.cpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_128_bgq.cpp \
	dirac_operators/stD/dirac_operator_stD_bgq.cpp
include_HEADERS+= \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_bgq.hpp \
	dirac_operators/tmQ/dirac_operator_tmQ_bgq.hpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_bgq.hpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_128_bgq.hpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_bgq.hpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ_bgq.hpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_bgq.hpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_128_bgq.hpp \
	dirac_operators/stD/dirac_operator_stD_bgq.hpp
endif

EXTRA_DIST= \
	base/git_info.hpp \
	dirac_operators/WclovQ/dirac_operator_WclovQ_portable.cpp \
	dirac_operators/stD/dirac_operator_stD_portable.cpp \
	dirac_operators/stD/dirac_operator_stD_32_portable.cpp \
	dirac_operators/stD/dirac_operator_stD_bgq_template.cpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_portable.cpp \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_portable.cpp \
	dirac_operators/tmQ/dirac_operator_tmQ_portable.cpp \
	dirac_operators/tmQ_left/dirac_operator_tmQ_left_portable.cpp
