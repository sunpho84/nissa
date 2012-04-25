GCC=gcc
CC=openmpicxx
CFLAGS=-O2 -Wall -DSVN_VERS=0
INCLUDE_PATH=src /Users/francesco/Prace/Programs/lemon/include /opt/local/include/openmpi/
LIBRARY_PATH=/Users/francesco/Prace/Programs/lemon/lib

VPATH=projects/semileptonic

################################################ define the programs ############################################

programs=projects/semileptonic/semileptonic_smeared projects/eight_BK/smeared_BK_all_in_one

################################################## global targets ###############################################

all: src/libnissa.a $(programs)

clean:
	 rm -rf $(addsuffix .d,$(nissa_library_pieces)) $(addsuffix .o,$(nissa_library_pieces)) src/libnissa.a
.PHONY: clean

############################################# define the library pieces #########################################

base=$(addprefix base/, close communicate debug global_variables init random routines vectors)
dirac_operators=$(addprefix dirac_operators/dirac_operator_, stD/dirac_operator_stD tmDeoimpr/dirac_operator_tmDeoimpr tmQ/dirac_operator_tmQ tmQ_left/dirac_operator_tmQ_left tmQ2/dirac_operator_tmQ2 tmQ/dirac_operator_tmQ_128 tmQ/reconstruct_tm_doublet tmQ2/dirac_operator_tmQ2_128)
geometry=$(addprefix geometry/, geometry_eo geometry_lx geometry_mix)
inverters=$(addprefix inverters/twisted_mass/, cg_invert_tmDeoimpr cg_invert_tmQ2 cg_128_invert_tmQ2 cgm_invert_tmQ2 tm_frontends) 
IO=$(addprefix IO/, checksum endianess input reader writer)
new_types=$(addprefix new_types/, complex dirac float128 rat_exp spin su3)
operations=$(addprefix operations/, contract fft fourier_transform gauge_fixing gaugeconf remap_vector smear su3_paths vector_gather)

#include all pieces in the library, define objects and target libary
nissa_library_pieces=$(addprefix src/, $(base) $(dirac_operators) $(inverters) $(geometry) $(IO) $(new_types) $(operations) linalgs/linalgs)
nissa_library_objects: $(addsuffix .o,$(nissa_library_pieces))

############################## rules to produce objects, library and programs ###################################

$(programs): src/libnissa.a Makefile
	$(CC) src/libnissa.a $(addprefix -I,$(INCLUDE_PATH)) $(addprefix -L,$(LIBRARY_PATH)) $(CFLAGS) -llemon -o $@ $(addsuffix .cpp, $@)

#$(addsuffix .d,$(programs)): %.d: %.cpp Makefile
#	$(GCC) $< -MM $(addprefix -I,$(INCLUDE_PATH)) -o $@

$(addsuffix .o,$(nissa_library_pieces)): %.o: %.cpp
	$(CC) $(CFLAGS) $< $(addprefix -I,$(INCLUDE_PATH)) -c -o $@

src/libnissa.a: nissa_library_objects
	ar cru src/libnissa.a $(addsuffix .o,$(nissa_library_pieces))

#dependencies:  $(addsuffix .d,$(programs))
#	@ echo "generate tables of dependencies"

