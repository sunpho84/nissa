GCC=gcc
CC=openmpicc
CFLAGS=-O4 -Wall -DSVN_VERS=0
INCLUDE_PATH=src ~/Prace/Programs/lemon/include 

programs=semileptonic_smeared

all: dependencies

#$(programs) : $(addsuffix .cpp,$(programs))
#	$(CC) $(INCLUDE)-c $^ $(CFLAGS) -o $@

$(addsuffix .d,$(programs)): %.d: %.cpp Makefile
	@ $(GCC) $< -MM $(addprefix -I,$(INCLUDE_PATH)) -o $@


dependencies:  $(addsuffix .d,$(programs))
	@ echo "generate tables of dependencies"

VPATH=projects/semileptonic

