#path to folder including nissa.h and libnissa.a
NISSA_PATH=
#name of the c++ compiler, e.g. mpicxx
CC=
#path where "lib" and "include" folders are present
LEMON_PATH=
#to be fit according to nissa Makefile
LEMON_TYPE=
#additional flags like -O2
CFLAGS=
#this can be tipicaly left to gcc
GCC=gcc

#derived parameters
SVN_VERS:=$(shell svnversion 2>/dev/null)
MACROS=$(addprefix -D,SVN_VERS=\"$(SVN_VERS)\" LEMON_TYPE=$(LEMON_TYPE))
INCLUDE_PATHS=$(LEMON_PATH)/include $(NISSA_PATH)
INCLUDE_FLAGS=$(addprefix -I,$(INCLUDE_PATHS))
LIBRARY_PATHS=$(LEMON_PATH)/lib $(NISSA_PATH)
LIBRARY_FLAGS=$(addprefix -L,$(LIBRARY_PATHS))

################################################## global targets ###############################################

all: semileptonic_smeared_point_source_tens

semileptonic_smeared_point_source_tens: semileptonic_smeared_point_source_tens.o $(NISSA_PATH)/libnissa.a
	$(CC) -o semileptonic_smeared_point_source_tens semileptonic_smeared_point_source_tens.o $(LIBRARY_FLAGS) -llemon -lnissa

semileptonic_smeared_point_source_tens.o: semileptonic_smeared_point_source_tens.cpp
	$(CC) $(CFLAGS) $(INCLUDE_FLAGS) -c -o semileptonic_smeared_point_source_tens.o semileptonic_smeared_point_source_tens.cpp

clean:
	rm -fv semileptonic_smeared_point_source_tens.o semileptonic_smeared_point_source_tens

.PHONY: clean