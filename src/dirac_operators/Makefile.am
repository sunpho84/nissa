AM_CPPFLAGS=-I$(top_srcdir)/src

#compile or not bgq test and source
if BGQ
bgq_src_files= \
        dirac_operator_tmQ/dirac_operator_tmQ_bgq.cpp \
        dirac_operator_tmQ2/dirac_operator_tmQ2_bgq.cpp \
        dirac_operator_tmQ/dirac_operator_tmQ_bgq.hpp \
        dirac_operator_tmQ2/dirac_operator_tmQ2_bgq.hpp \
	dirac_operator_tmQ2/dirac_operator_tmQ2_128_bgq.cpp \
	dirac_operator_tmQ2/dirac_operator_tmQ2_128_bgq.hpp \
	dirac_operator_stD/dirac_operator_stD_bgq.cpp \
	dirac_operator_stD/dirac_operator_stD_bgq.hpp
endif


noinst_LTLIBRARIES=libdirac_operators.la
libdirac_operators_la_SOURCES= \
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

dirac_operators_stD_includedir=$(includedir)/dirac_operators/dirac_operator_stD
dirac_operators_stD_include_HEADERS=dirac_operator_stD/dirac_operator_stD.hpp

dirac_operators_tmDeo_impr_includedir=$(includedir)/dirac_operators/dirac_operator_tmDeoimpr
dirac_operators_tmDeo_impr_include_HEADERS= \
	dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr.hpp \
        dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr_128.hpp

dirac_operators_tmQ_includedir=$(includedir)/dirac_operators/dirac_operator_tmQ
dirac_operators_tmQ_include_HEADERS= \
	dirac_operator_tmQ/dirac_operator_tmQ.hpp \
        dirac_operator_tmQ/reconstruct_tm_doublet.hpp \
        dirac_operator_tmQ/dirac_operator_tmQ_128.hpp

dirac_operators_tmQ_left_includedir=$(includedir)/dirac_operators/dirac_operator_tmQ_left
dirac_operators_tmQ_left_include_HEADERS=dirac_operator_tmQ_left/dirac_operator_tmQ_left.hpp

dirac_operators_tmQ2_includedir=$(includedir)/dirac_operators/dirac_operator_tmQ2
dirac_operators_tmQ2_include_HEADERS=dirac_operator_tmQ2/dirac_operator_tmQ2.hpp dirac_operator_tmQ2/dirac_operator_tmQ2_128.hpp

dirac_operators_tmclovQ_includedir=$(includedir)/dirac_operators/dirac_operator_tmclovQ
dirac_operators_tmclovQ_include_HEADERS= \
	dirac_operator_tmclovQ/dirac_operator_tmclovQ.hpp \
	dirac_operator_tmclovQ/reconstruct_tmclov_doublet.hpp

dirac_operators_Wstat_includedir=$(includedir)/dirac_operators/dirac_operator_Wstat
dirac_operators_Wstat_include_HEADERS=dirac_operator_Wstat/dirac_operator_Wstat.hpp

dirac_operators_WclovQ_includedir=$(includedir)/dirac_operators/dirac_operator_WclovQ
dirac_operators_WclovQ_include_HEADERS=dirac_operator_WclovQ/dirac_operator_WclovQ.hpp

dirac_operators_WclovQ2_includedir=$(includedir)/dirac_operators/dirac_operator_WclovQ2
dirac_operators_WclovQ2_include_HEADERS=dirac_operator_WclovQ2/dirac_operator_WclovQ2.hpp

dirac_operators_tmclovQ2_includedir=$(includedir)/dirac_operators/dirac_operator_tmclovQ2
dirac_operators_tmclovQ2_include_HEADERS=dirac_operator_tmclovQ2/dirac_operator_tmclovQ2.hpp

EXTRA_DIST= \
        dirac_operator_WclovQ/dirac_operator_WclovQ_portable.cpp \
        dirac_operator_stD/dirac_operator_stD_portable.cpp \
        dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr_portable.cpp \
	dirac_operator_tmQ/dirac_operator_tmQ_portable.cpp \
        dirac_operator_tmQ_left/dirac_operator_tmQ_left_portable.cpp
