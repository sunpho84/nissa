#externally defined path and variables
CC=openmpicxx
LEMON_PATH=/Users/francesco/Prace/Programs/lemon
LEMON_TYPE=uint64_t
CFLAGS=-ggdb3

#derived parameters
SVN_VERS:=$(shell svnversion 2>/dev/null)
MACROS=$(addprefix -D,SVN_VERS=\"$(SVN_VERS)\" LEMON_TYPE=$(LEMON_TYPE))
INCLUDE_PATH=../../src $(addsuffix /include,$(LEMON_PATH))
LIBRARY_PATH=$(addsuffix /lib,$(LEMON_PATH))

GCC=gcc

################################################ define the objects ############################################

#include all objects
objects=twisted_propagator

############################################# define the executables ############################################

#collect all the executables
executables=test_twisted_propagator

################################################## global targets ###############################################

all: Makefile $(executables)

clean:
	 rm -rf $(addsuffix .d,$(objects) $(executables)) $(addsuffix .o,$(objects)) $(executables)

#include table of dependencies
-include $(addsuffix .d,$(projects) $(executables))

############################################ rules to produce objects ############################################

$(addsuffix .o,$(objects)): %.o: %.cpp %.d Makefile
	$(CC)				\
	$(addprefix -I,$(INCLUDE_PATH)) \
	$<				\
	$(CFLAGS) 			\
	$(MACROS) 			\
	-c -o $@

######################################### rules to produce dependencies ##########################################

$(addsuffix .d,$(projects) $(executables)): %.d: %.cpp Makefile
	$(GCC) $< -MM -MT $(@:%.d=%.o) $(addprefix -I,$(INCLUDE_PATH)) -o $@

######################################### rules to produce executables ############################################

$(executables): %: %.cpp %.d ../../src/libnissa.a Makefile
	$(CC)                           \
	$(addprefix -I,$(INCLUDE_PATH)) \
	$(addprefix -L,$(LIBRARY_PATH)) \
	$(CFLAGS)                       \
	$(MACROS)                       \
	-o $@                           \
	$(addsuffix .cpp, $@)           \
	$(addsuffix .o, $(objects))     \
	../../src/libnissa.a		\
	-llemon


############################################## phonyfyse the needed target ##########################################

.PHONY: clean
