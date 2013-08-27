#compile or not bgq test and source
if BGQ
bgq_src_files= \
        dirac_operator_tmQ/dirac_operator_tmQ_bgq.cpp \
        dirac_operator_tmQ/dirac_operator_tmQ_bgq.h \
        dirac_operator_tmQ2/dirac_operator_tmQ2_bgq.cpp \
        dirac_operator_tmQ2/dirac_operator_tmQ2_bgq.h
else
bgq_src_files=
endif


noinst_LTLIBRARIES=libdirac_operators.la
libdirac_operators_la_SOURCES= \
        dirac_operator_stD/dirac_operator_stD.h \
        dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr.h \
        dirac_operator_tmQ/dirac_operator_tmQ.h \
        dirac_operator_tmQ_left/dirac_operator_tmQ_left.h \
        dirac_operator_tmQ2/dirac_operator_tmQ2.h \
        dirac_operator_tmQ/reconstruct_tm_doublet.h \
        dirac_operator_tmclovQ/reconstruct_tmclov_doublet.h \
        dirac_operator_Wstat/dirac_operator_Wstat.h \
        dirac_operator_WclovQ/dirac_operator_WclovQ.h \
        dirac_operator_WclovQ2/dirac_operator_WclovQ2.h \
        dirac_operator_tmclovQ/dirac_operator_tmclovQ.h \
        dirac_operator_tmclovQ2/dirac_operator_tmclovQ2.h \
        dirac_operator_tmQ/dirac_operator_tmQ_128.h \
        dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr_128.h \
        dirac_operator_tmQ2/dirac_operator_tmQ2_128.h \
        dirac_operator_stD/dirac_operator_stD.cpp \
        dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr.cpp \
        dirac_operator_tmQ/dirac_operator_tmQ.cpp \
        dirac_operator_tmQ_left/dirac_operator_tmQ_left.cpp \
        dirac_operator_tmQ2/dirac_operator_tmQ2.cpp \
        dirac_operator_tmQ/reconstruct_tm_doublet.cpp \
        dirac_operator_tmclovQ/reconstruct_tmclov_doublet.cpp \
        dirac_operator_Wstat/dirac_operator_Wstat.cpp \
        dirac_operator_WclovQ/dirac_operator_WclovQ.cpp \
        dirac_operator_WclovQ2/dirac_operator_WclovQ2.cpp \
        dirac_operator_tmclovQ/dirac_operator_tmclovQ.cpp \
        dirac_operator_tmclovQ2/dirac_operator_tmclovQ2.cpp \
        dirac_operator_tmQ/dirac_operator_tmQ_128.cpp \
        dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr_128.cpp \
        dirac_operator_tmQ2/dirac_operator_tmQ2_128.cpp \
	$(bgq_src_files)