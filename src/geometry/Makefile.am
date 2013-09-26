AM_CPPFLAGS=-I$(top_srcdir)/src

#compile or not virtual parallelization
if USE_VNODES
vnodes_src_files=geometry_vir.cpp
vnodes_header_files=geometry_vir.hpp
endif

noinst_LTLIBRARIES=libgeometry.la
libgeometry_la_SOURCES=$(vnodes_src_files) geometry_eo.cpp geometry_lx.cpp geometry_mix.cpp geometry_Wsklx.cpp

geometry_includedir=$(includedir)/geometry
geometry_include_HEADERS=$(vnodes_header_files) geometry_eo.hpp geometry_lx.hpp geometry_mix.hpp geometry_Wsklx.hpp
