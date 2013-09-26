AM_CPPFLAGS=-I$(top_srcdir)/src

noinst_LTLIBRARIES=libhmc.la
libhmc_la_SOURCES= \
        backfield.cpp \
        gauge/gluonic_force.cpp \
        gauge/tree_level_Symanzik_force.cpp \
        gauge/tree_level_Symanzik_action.cpp \
        gauge/Wilson_force.cpp \
        momenta/momenta_generation.cpp \
        momenta/momenta_action.cpp \
        rootst_eoimpr/adaptative_control.cpp \
        rootst_eoimpr/rootst_eoimpr_action.cpp \
        rootst_eoimpr/rootst_eoimpr_eigenvalues.cpp \
        rootst_eoimpr/rootst_eoimpr_quark_force.cpp \
        rootst_eoimpr/rootst_eoimpr_omelyan_integrator.cpp \
        rootst_eoimpr/rootst_eoimpr_pseudofermions_generation.cpp \
        rootst_eoimpr/rootst_eoimpr_rhmc_step.cpp

hmc_includedir=$(includedir)/hmc
hmc_include_gaugedir=$(includedir)/hmc/gauge
hmc_include_momentadir=$(includedir)/hmc/momenta
hmc_include_rootst_eoimprdir=$(includedir)/hmc/rootst_eoimpr

hmc_include_HEADERS= backfield.hpp
hmc_include_gauge_HEADERS= \
        gauge/gluonic_force.hpp \
        gauge/tree_level_Symanzik_force.hpp \
        gauge/tree_level_Symanzik_action.hpp \
        gauge/Wilson_force.hpp
hmc_include_momenta_HEADERS= \
        momenta/momenta_generation.hpp \
        momenta/momenta_action.hpp
hmc_include_rootst_eoimpr_HEADERS= \
        rootst_eoimpr/adaptative_control.hpp \
        rootst_eoimpr/rootst_eoimpr_action.hpp \
        rootst_eoimpr/rootst_eoimpr_eigenvalues.hpp \
        rootst_eoimpr/rootst_eoimpr_quark_force.hpp \
        rootst_eoimpr/rootst_eoimpr_omelyan_integrator.hpp \
        rootst_eoimpr/rootst_eoimpr_pseudofermions_generation.hpp \
        rootst_eoimpr/rootst_eoimpr_rhmc_step.hpp

EXTRA_DIST= rootst_eoimpr/rat_expansion_database.cpp
