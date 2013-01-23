tmLQCD Optimization for BlueGene/Q

Speed
=====

Showcase performance per node (204000 machine peak, 159120 theoretical limit, 131000 practical limit):

[tmLQCD_official_373, tmLQCD_official_375]
64042 MFlop/s (12x10x10x12 local lattice, 64 threads, 4x4x2x1=32 nodes, double precision, 1608 Flop/stencil) 
85840 MFlop/s (12x10x10x12 local lattice, 64 threads, 4x4x2x1=32 nodes, single precision, 1608 Flop/stencil)
100287 MFlop/s (12x10x10x12 local lattice, 64 threads, without communication, double precision, 1608 Flop/stencil)
110335 MFlop/s (12x10x10x12 local lattice, 64 threads, without communication, single precision, 1608 Flop/stencil)
53261 MFlop/s (12x10x10x12 local lattice, 64 threads, 4x4x2x1=32 nodes, single precision, 1320 Flop/stencil)
87119 MFlop/s (12x10x10x12 local lattice, 64 threads, without communication, double precision, 1320 Flop/stencil)
99642 MFlop/s (12x10x10x12 local lattice, 64 threads, without communication, single precision, 1320 Flop/stencil)

75978 MFlop/s (12x10x10x12 local lattice, 48 threads,  4x4x4x1=64 nodes, double precision, 1608 Flop/stencil)
90163 MFlop/s (12x10x10x12 local lattice, 64 threads,  4x4x4x1=64 nodes, single precision, 1608 Flop/stencil)

[tmLQCD_official_450]
60113 MFlop/s (12x10x10x12 local lattice, 64 threads, 4x4x2x1=32 nodes, benchmark.c)
114404 MFlop/s  (12x10x10x12 local lattice, 64 threads, 4x4x2x1=32 nodes, benchmark.c, Communication switched off)

[tmLQCD_official_415]
50387 MFlop/s (12x10x8x8 local lattice, 64 threads, 4x4x2x1=32 nodes, double precision, 1608 Flop/stencil)
94634 MFlop/s (12x10x8x8 local lattice, 64 threads, no communication, double precision, 1608 Flop/stencil)

[tmLQCD_official_378]
46364 MFlop/s (16x8x4x8 local lattice, 64 threads, 4x4x4x1=64 nodes, double precision, 1608 Flop/stencil)
84608 MFlop/s (16x8x4x8 local lattice, 64 threads, no communication, double precision, 1608 Flop/stencil)

[tmLQCD_official_379]
18175 MFlop/s (8x16x2x4 local lattice, 64 threads, 4x4x4x1=64 nodes, double precision, 1608 Flop/stencil)
60758 MFlop/s (8x16x2x4 local lattice, 64 threads, no communication, double precision, 1608 Flop/stencil)

[tmLQCD_official_465]
50440 MFlop/s (12x8x8x8 local lattice, 64 threads, 4x4x2x1=32 nodes, double precision, 1608 Flop/stencil)
92730 MFlop/s (12x8x8x8 local lattice, 64 threads, no communication, double precision, 1608 Flop/stencil)

[tmLQCD_official_469]
23286 MFlop/s (24x4x4x4 local lattice, 64 threads, 4x4x2x1=32 nodes, double precision, 1608 Flop/stencil)
53938 MFlop/s (24x4x4x4 local lattice, 64 threads, no communication, double precision, 1608 Flop/stencil)

[tmLQCD_official_474]
65293 MFlop/s (12x10x10x12 local lattice, 64 threads, 4x4x4x1=64 nodes, EDCBAT mapping, double precision, 1608 Flop/stencil)

[tmLQCD_official_476]
75359 MFlop/s (12x10x10x12 local lattice, 48 threads, 4x4x4x1=64 nodes, double precision, 1608 Flop/stencil)
87737 MFlop/s (12x10x10x12 local lattice, 64 threads, 4x4x4x1=64 nodes, single precision, 1608 Flop/stencil)


Best performance is reached with square big local lattices, but whose working set still fits into the L2 cache. This limit is at about 12x10x10x12. For bigger local lattices, kernel performance drops by 74%.


Compile
=======

I did not adapt the build system, just fixed some issues. To compile, there are two scripts: build-bgq.sh used to compile on JUQUEEN and build-x86.sh for my desktop computer I use for debugging.
Descriptions of additional preprocessor flags:
NDEBUG			Defining it removes debugging instructions including asserts
NVALGRIND		Define if valgrind.h is not available
BGQ_QPX			BGQ_QPX=1 to use BG/Q unique features like AXU instructions, BGQ_QPX=0 to emulate them
PAPI			To enable BGPM hardware performance counter diagnostics in bgqbench program (PAPI is actually not used, has been in a previous version)
BGQ_UNVECTORIZE	BGQ_UNVECTORIZE=1 will extract the relevant data from qpx vectors before sending it in T-direction; BGQ_UNVECTORIZE=1 just sends the complete vector and therefore sends twice as much data in T-direction, but saves an intermediate step.
BGQ_COORDCHECK	Define to store the expected coordinate into the fields instead the dloation point data. Used to check correct indexing.
BGQ_REPLACE		Replace many operators by its BG/Q-optimized equivialent

./configure options:
--enable-gaugecopy			Gaugecopy is always used. Do not disable.
--with-mpidimension=XYZT	Always use 4 dimension, actual dimensions is determined at runtime depending on the topology
--enable-optimize=no		There is currently no distinction betweeb MODULES and SMODULES
--enable-spi				Speaks for itself; MPI will be used if disabled
--disable-omp				OpenMP and the BGQ mechanism interfere (nested parallelism), although it somehow works.

Some notable build system issues were:
- Executables are always linked which may take a lot of time if Link-Time Optimization is enabled. The issue was that some .PHONY targets are in the dependency lists.
- GIT-VERSION-GEN has not been called when something changes. If called, it's called eveytime if the git version did not change because the timestamp of git_hash.h is not updated
- Non-consistent use of CPPFLAGS, CFLAGS, DEPFLAGS, LDFLAGS, LIBS
- AC_FUNC_MALLOC may change all calls to malloc to calls to rpl_malloc, which must be defined somewhere


Implementation Notes
====================

Layouts
-------

The new spinorfield type is bgq_weylfield_controlblock. As its name says, it does not contain the data itself but information how to find it. Data can be one of 5 layouts, their memory is allocated on demand:
ly_full_double	Default layout, 192 byte per spinor.
ly_full_float	Same as above, but in single precision, 96 byte per spinor
ly_weyl_double	Halfspinor layout (one halfspinor per direction) as written by HoppingMatrix, 768 byte per spinor; The reading function has to compute the spinor
ly_weyl_float	Same as above, but in single precision, 384 byte per spinor
ly_legacy		For compatibility with existing functions, double precision, non-vectorized

The ordering of spinors in ly_legacy is "eosub", for all others it is "collapsed", where the spinors on the communication surface come first and two spinors in t-direction are combined to one physical site to allow optimial vectorization.

A specific layout can be requested using the function bgq_spinorfield_prepareWrite. The function bgq_spinorfield_prepareRead also translates the existing data from some other layout. Because the order depends on whether the field contains odd or even spinors, this information must be known for every field operation.

"even_odd_flag" must be defined. Non-even/odd spinor fields are not supported.

For HoppingMatrix, there another arrangement of the gauge field is required, which is handled like _GAUGE_COPY. It is always kept in double precesion, even if single precision is chosen.


Local Lattice Constraints
-------------------------

The chosen memory layout puts some constraints on the size of local lattice (lengths T, LX, LY and LZ). Or have to be special-cased. These are:

- The size of every communication buffer must be a multiple of 4 (i.e. T*LX*LY, T*LX*LZ, T*LY*LZ and LX*LY*LZ in case of 4-dimensional MPI/SPI-parallelization)
  (This is due to SPI requirering buffers to be 128-byte-aligned. This can be special-cased by adding some padding)

- LX,LY,LZ must each be multiples of 2
  (Violating this means that local lattices on different nodes have different even/odd-geometry, which is not implemented)

- T must be a multiple 4
  (One physical site stores either two even or two odd sites, therefore logically spanning 4 sites)

- T should be at least 8
  (T=4 means that a physical site neighbours both local borders, while I see no reason it should not work, I did not try it)

- (T-4)*(LX-2)*(LY-2)*(LZ-2) / T*LX*LY*LZ, i.e. the share of the local value that is not at the surface, should be as big as possible. In case of less than 4-dimensional distributed parallelization, remove the minuend from the non-parallelized dimension.
  (The body can be computed while the surface is transfering. The larger the body, the more computation can be hidden behind communication)

Take care, these constraints are not checked in the program.


Symmetric Multithreading
------------------------

OpenMP is used to spawn threads. Everything else is done using a custom mechanism.
The function bgq_parallel starts multithreading and calls a function that runs only on one thread.
This master thread then may call bgq_master_call to call a worker function using all threads, including master.
bgq_master_sync called in the master thread waits for all worker threads to finish.

The overhead per bgq_master_call followed by bgq_master_sync is about 3000 cycles for 64 threads.

bgq_parallel_mainlike can be used to use a main(int argc, char *argv[]) as the controlling master thread.

The typical pattern for workers is to accept a pointer to a static struct containing the parameters, including the amount of work to do. Using the thread-id (tid) and the total number of available threads, the worker computes its part of the total workload. The master should, before setting up the parameter struct, call bgq_master_sync because there might be a worker from the previous round still accessing it. 


Reductions and Spinor-wise Operators
-------------------------------------

To avoid an extra-translation step, operators should read the data in whatever layout is available, without using conditionals in the kernel iteration. There are 5 layouts, therefore binary operators have 50 combinations (25 times 2 for writing in single/double precision). The files bgq_reduction.inc.c and bgq_operator.inc.c contain the boilerplate code for this. It relies on the inlining capitibility of the compiler. To avoid excessive compilation times, not all combinations actually exist.

bgq_stdreductions.c and bgq_stdoperators.c contain some common operators.

Important note: bgq_operator.inc.c contains inline assembler to write the result spinor without memory clobber or reorder barrier. If input- and aoutput-field are identical, the compiler might reorder the writes before reading it. However, for the predefined operators, there is a data dependency beween them.


Legacy Operators
----------------

Not all field functions have been translated to use the new layout. For the remaining, bgq_translate_spinorfield can be used to get the bgq_weylfield_controlblock that corresponds to a legacy field. Use bgq_weylfield_controlblock.legacy_field for the other direction.

To be sure that the legacy field is up-to-date, spinorfield_enable must be called before accessing it. This will copy the field data into the expected location, therefore involves a performance penalty.

Field operatores that have been disabled using #if !BGQ_REPLACE, and their replacement resides in bgq/bgq_legacy.c . This effectively leaves the files empty and the compiler might warn about it.

Notably, deriv_Sb has not been replaced.


Error-Checking
--------------

assert() is the primary statement used for checking errornous conditions. Try undefining NDEBUG preprocessor macro if something does not seem right. Non-performance relevant asserts might be converted to permanent error-checking conditionals.


Freeing Resources
-----------------

None of the resource allocated are ever freed. This is due to the nature of the application possibly needing the resources until it finishes, i.e. when the OS will free everything anyways.

However some fields, noteably those used by solvers, are allocated and freed at when they are no more needed. These fields are put into a pool to be reused when fields are allocated again. This saves runtime because pointer locations need to be set up for HoppingMatrix.


Asynchronous Communication
--------------------------

Communication overlaps with the processing of sites that do not need any data from neighbor nodes. In addition, communication is not required to finish before
- Field field result is used in some field operator OR
- Another communication is triggered that needs the communication buffers

Therefore communication can be hidded even more if, after HoppingMatrix has been called, a linalg operator on a different field should be called. E.g. 
HoppingMatrix(EO, targetfield, source)
add(temp, source, anothersource)
dosomething(targetfield)


More Optimization Possibilities
-------------------------------

If the working set does not fit into the L2 cache, a full spinor optimization will be faster.

The prefetching depth (dcbt) has not been fine-tuned.

Currently tm_times_Hopping_Matrix is split up into to operations, using in a single loop is probably faster. More generally, any field could be annotated with a factor to multiply with on the fly. This, however, will again increase the number of combinations of fields layouts the operators must implement. tm_sub_Hopping_Matrix is probably not worth doing this because it reads another stream which the hardware is unable to track (i.e. 4 prefetch streams or more).


bgqbench
--------

There is a new benchmarking program bgqbench.c especially designed to compare the execution speed of HoppingMatrix. It shows the speed with various number of threads and options, like double/single precision, SPI/MPI/no communication and prefetching modes.

Set BGQ_REPLACE=0 such that it can show the error compared to the legacy implementation.

