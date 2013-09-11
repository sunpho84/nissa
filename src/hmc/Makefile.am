AM_CXXFLAGS=-I$(top_srcdir)/src

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

hmcincludedir=$(includedir)/hmc
hmcinclude_HEADERS= \
	backfield.h \
        gauge/gluonic_force.h \
        gauge/tree_level_Symanzik_force.h \
        gauge/tree_level_Symanzik_action.h \
        gauge/Wilson_force.h \
        momenta/momenta_generation.h \
        momenta/momenta_action.h \
        rootst_eoimpr/adaptative_control.h \
        rootst_eoimpr/rootst_eoimpr_action.h \
        rootst_eoimpr/rootst_eoimpr_eigenvalues.h \
        rootst_eoimpr/rootst_eoimpr_quark_force.h \
        rootst_eoimpr/rootst_eoimpr_omelyan_integrator.h \
        rootst_eoimpr/rootst_eoimpr_pseudofermions_generation.h \
        rootst_eoimpr/rootst_eoimpr_rhmc_step.h

EXTRA_DIST= \
	rootst_eoimpr/rat_expansion_database.cpp