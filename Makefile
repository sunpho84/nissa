GCC=gcc
CC=openmpicc
CFLAGS=-O2 -Wall -DSVN_VERS=0
INCLUDE_PATH=src /Users/francesco/Prace/Programs/lemon/include /opt/local/include/openmpi/

new_types=$(addprefix new_types/, complex dirac float128 rat_exp spin su3)
modules=$(addprefix src/,$(new_types))

all: objects

#$(programs) : $(addsuffix .cpp,$(programs))
#	$(CC) $(INCLUDE)-c $^ $(CFLAGS) -o $@

$(addsuffix .d,$(programs)): %.d: %.cpp Makefile
	@ $(GCC) $< -MM $(addprefix -I,$(INCLUDE_PATH)) -o $@

$(addsuffix .o,$(modules)): %.o: %.cpp Makefile
	@ $(CC) $< $(addprefix -I,$(INCLUDE_PATH)) -c -o $@


#dependencies:  $(addsuffix .d,$(programs))
#	@ echo "generate tables of dependencies"

objects: $(addsuffix .o,$(modules))

VPATH=projects/semileptonic

clean:
	@ rm -rvf $(addsuffix .d,$(modules)) $(addsuffix .o,$(modules))
.PHONY: clean
