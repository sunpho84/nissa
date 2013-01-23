AUTOMAKE_OPTIONS=subdir-objects

#add where to find the library and tell to link
AM_LDFLAGS=-L@top_builddir@ -L/Users/francesco/Prace/nissa/
LDADD=@LIBS@ libinterface.a ../../libnissa.a

#target library nissa
lib_LIBRARIES=libinterface.a libkruse.a

####################################################################################################################################################

libinterface_a_SOURCES=src/interface/additional_routines.c src/interface/additional_variables.c 
libkruse_a_SOURCES=src/kruse/bgq_HoppingMatrix.c src/kruse/bgq_gaugefield.c

#list of all the programs
bin_PROGRAMS=main

#program to generate confs
main_SOURCES=pseudo_main.cpp
