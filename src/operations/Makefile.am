AM_CXXFLAGS=-I$(top_srcdir)/src

noinst_LTLIBRARIES=liboperations.la
liboperations_la_SOURCES= \
        contract/mesons_2pts.cpp \
        contract/mesons_eight.cpp \
        contract/optimized_mesons_2pts.cpp \
        contract/site_contract.cpp \
        contract/stag.cpp \
        covariant_derivative.cpp \
        fft.cpp \
        fourier_transform.cpp \
        gauge_fixing.cpp \
        gaugeconf.cpp \
        remap_vector.cpp \
        remez/remez_algorithm.cpp \
        shift.cpp \
        source.cpp \
        smearing/APE.cpp \
        smearing/HYP.cpp \
        smearing/stout.cpp \
        smearing/gaussian.cpp \
        su3_paths/all_rectangles.cpp \
        su3_paths/arbitrary.cpp \
	su3_paths/hypercubic_transfer.cpp \
        su3_paths/plaquette.cpp \
        su3_paths/pline.cpp \
        su3_paths/rectangles.cpp \
        su3_paths/rectangular_staples.cpp \
        su3_paths/shortest_hypercubic_paths.cpp \
        su3_paths/squared_staples.cpp \
        su3_paths/topological_charge.cpp \
        vector_gather.cpp

operations_includedir=$(includedir)/operations
operations_include_HEADERS= \
        covariant_derivative.h \
        fft.h \
        fourier_transform.h \
        gauge_fixing.h \
        gaugeconf.h \
        remap_vector.h \
        shift.h \
        source.h \
        vector_gather.h

operations_remez_includedir=$(includedir)/operations/remez
operations_remez_include_HEADERS=remez/remez_algorithm.h

operations_contract_includedir=$(includedir)/operations/contract
operations_contract_include_HEADERS= \
        contract/mesons_2pts.h \
        contract/mesons_eight.h \
        contract/optimized_mesons_2pts.h \
        contract/site_contract.h \
        contract/stag.h

operations_smearing_includedir=$(includedir)/operations/smearing
operations_smearing_include_HEADERS= \
        smearing/APE.h \
        smearing/HYP.h \
        smearing/stout.h \
        smearing/gaussian.h

operations_su3_paths_includedir=$(includedir)/operations/su3_paths
operations_su3_paths_include_HEADERS= \
        su3_paths/all_rectangles.h \
        su3_paths/arbitrary.h \
	su3_paths/hypercubic_transfer.h \
        su3_paths/plaquette.h \
        su3_paths/pline.h \
        su3_paths/rectangles.h \
        su3_paths/rectangular_staples.h \
        su3_paths/shortest_hypercubic_paths.h \
        su3_paths/squared_staples.h \
        su3_paths/topological_charge.h
