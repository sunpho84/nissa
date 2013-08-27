#compile or not virtual parallelization
if USE_VNODES
vnodes_src_file=geometry_vir.h geometry_vir.cpp
endif
eo_src_file=geometry_eo.h geometry_eo.cpp
lx_src_file=geometry_lx.h geometry_lx.cpp
mix_src_file=geometry_mix.h geometry_mix.cpp
Wsklx_src_file=geometry_Wsklx.h geometry_Wsklx.cpp

noinst_LTLIBRARIES=libgeometry.la
libgeometry_la_SOURCES=$(eo_src_file) $(lx_src_file) $(mix_src_file) $(Wsklx_src_file) $(vnodes_src_file)