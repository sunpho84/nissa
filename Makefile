GCC=gcc
CC=openmpicc
CFLAGS=-O2 -Wall -DSVN_VERS=0
INCLUDE_PATH=src /Users/francesco/Prace/Programs/lemon/include /opt/local/include/openmpi/

#define the library pieces
geometry=$(addprefix geometry/, geometry_eo geometry_lx geometry_mix)
IO=$(addprefix IO/, checksum endianess reader writer)
new_types=$(addprefix new_types/, complex dirac float128 rat_exp spin su3)
operations=$(addprefix operations/, contract fft fourier_transform gauge_fixing gaugeconf remap_vector smear su3_paths vector_gather)

modules=$(addprefix src/, $(geometry) $(IO) $(new_types) $(operations) linalgs/linalgs)

all: objects

#$(programs) : $(addsuffix .cpp,$(programs))
#	$(CC) $(INCLUDE)-c $^ $(CFLAGS) -o $@

$(addsuffix .d,$(programs)): %.d: %.cpp Makefile
	$(GCC) $< -MM $(addprefix -I,$(INCLUDE_PATH)) -o $@

$(addsuffix .o,$(modules)): %.o: %.cpp Makefile
	$(CC) $< $(addprefix -I,$(INCLUDE_PATH)) -c -o $@


#dependencies:  $(addsuffix .d,$(programs))
#	@ echo "generate tables of dependencies"

objects: $(addsuffix .o,$(modules))

VPATH=projects/semileptonic

clean:
	 rm -rf $(addsuffix .d,$(modules)) $(addsuffix .o,$(modules))
.PHONY: clean
