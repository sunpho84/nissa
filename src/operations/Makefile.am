AM_CPPFLAGS=-I$(top_srcdir)/src

#library name
noinst_LTLIBRARIES=liboperations.la

#sources
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

#headers
operations_includedir=$(includedir)/operations
operations_include_HEADERS= \
        covariant_derivative.hpp \
        fft.hpp \
        fourier_transform.hpp \
        gauge_fixing.hpp \
        gaugeconf.hpp \
        remap_vector.hpp \
        shift.hpp \
        source.hpp \
        vector_gather.hpp

operations_remez_includedir=$(includedir)/operations/remez
operations_remez_include_HEADERS=remez/remez_algorithm.hpp

operations_contract_includedir=$(includedir)/operations/contract
operations_contract_include_HEADERS= \
        contract/mesons_2pts.hpp \
        contract/mesons_eight.hpp \
        contract/optimized_mesons_2pts.hpp \
        contract/site_contract.hpp \
        contract/stag.hpp

operations_smearing_includedir=$(includedir)/operations/smearing
operations_smearing_include_HEADERS= \
        smearing/APE.hpp \
        smearing/HYP.hpp \
        smearing/stout.hpp \
        smearing/gaussian.hpp

operations_su3_paths_includedir=$(includedir)/operations/su3_paths
operations_su3_paths_include_HEADERS= \
        su3_paths/all_rectangles.hpp \
        su3_paths/arbitrary.hpp \
	su3_paths/hypercubic_transfer.hpp \
        su3_paths/plaquette.hpp \
        su3_paths/pline.hpp \
        su3_paths/rectangles.hpp \
        su3_paths/rectangular_staples.hpp \
        su3_paths/shortest_hypercubic_paths.hpp \
        su3_paths/squared_staples.hpp \
        su3_paths/topological_charge.hpp
