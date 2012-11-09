#externally defined path and variables
CC=openmpicxx
CFLAGS=-O0 -ggdb3

#derived parameters
SVN_VERS:=$(shell svnversion 2>/dev/null)
MACROS= #$(addprefix -D,SVN_VERS=\"$(SVN_VERS)\")
INCLUDE_PATH=../../src
LIBRARY_PATH=

GCC=gcc

################################################ define the objects ############################################

#include all objects
objects=

############################################# define the executables ############################################

#full programs
progs=Wstat_check point_Wstat_check

#collect all the executables
executables=$(tests) $(progs)

################################################## global targets ###############################################

all: Makefile $(executables)

clean:
	 rm -rf $(addsuffix .d,$(objects) $(executables)) $(addsuffix .o,$(objects)) $(executables)

#include table of dependencies
-include $(addsuffix .d,$(projects) $(executables))

############################################ rules to produce objects ############################################

$(addsuffix .o,$(objects) $(executables)): %.o: %.cpp %.d Makefile
	$(CC)				\
	$(addprefix -I,$(INCLUDE_PATH)) \
	$<				\
	$(CFLAGS) 			\
	$(MACROS) 			\
	-c -o $@

######################################### rules to produce dependencies ##########################################

$(addsuffix .d,$(objects) $(executables)): %.d: %.cpp Makefile
	$(GCC) $< -MM -MT $(@:%.d=%.o) $(addprefix -I,$(INCLUDE_PATH)) -o $@

######################################### rules to produce executables ############################################

$(executables): %: %.o $(addsuffix .o,$(objects)) ../../libnissa.a Makefile
	$(CC)                           \
	$(addprefix -I,$(INCLUDE_PATH)) \
	$(addprefix -L,$(LIBRARY_PATH)) \
	$(CFLAGS)                       \
	$(MACROS)                       \
	-o $@                           \
	$(addsuffix .o, $@)           \
	$(addsuffix .o, $(objects))     \
	../../libnissa.a

############################################## phonyfyse the needed target ##########################################

.PHONY: clean
