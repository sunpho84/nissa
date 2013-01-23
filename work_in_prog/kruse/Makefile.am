AUTOMAKE_OPTIONS=subdir-objects

#add where to find the library and tell to link
AM_LDFLAGS=-L@top_builddir@ -L/Users/francesco/Prace/nissa/
LDADD=@LIBS@ libinterface.a ../../libnissa.a

#target library nissa
lib_LIBRARIES=libinterface.a libkruse.a

####################################################################################################################################################

libinterface_a_SOURCES=src/interface/additional_routines.c src/interface/additional_variables.c 
libkruse_a_SOURCES=	\
	src/kruse/bgq_workers.inc.c \
	src/kruse/bgq_workers_double.c \
	src/kruse/bgq_workers_float.c \
	src/kruse/bgq_ReadWeyllayout.inc.c \
	src/kruse/bgq_HoppingMatrixWorker.inc.c \
	src/kruse/bgq_ComputeWeyl.inc.c \
	src/kruse/bgq_stdreductions.c \
	src/kruse/bgq_reduction.inc.c \
	src/kruse/bgq_stdoperators.c \
	src/kruse/bgq_operator.inc.c \
	src/kruse/asm1.c \
	src/kruse/asm2.c \
	src/kruse/bgq_qpx.c \
	src/kruse/bgq_gaugefield.c \
	src/kruse/bgq_spinorfield.c \
	src/kruse/bgq_comm.c \
	src/kruse/bgq_utils.c \
	src/kruse/bgq_field.c \
	src/kruse/bgq_HoppingMatrix.c \
	src/kruse/bgq_dispatch.c \
	src/kruse/mypapi.c

#bgq_legacy bgq_workers.inc bgq_workers_double bgq_workers_float bgq_ReadWeyllayout.inc bgq_HoppingMatrixWorker.inc bgq_ComputeWeyl.inc bgq_stdreductions bgq_reduction.inc bgq_stdoperators bgq_operator.inc asm1 asm2 bgq_qpx bgq_gaugefield bgq_spinorfield bgq_comm bgq_utils bgq_field bgq_HoppingMatrix bgq_dispatch mypapi


#list of all the programs
bin_PROGRAMS=main

#program to generate confs
main_SOURCES=pseudo_main.cpp
