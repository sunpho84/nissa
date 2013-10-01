AM_CPPFLAGS=-I$(top_srcdir)/src

#library name
noinst_LTLIBRARIES=libgeometry.la

#sources
libgeometry_la_SOURCES= \
	geometry_eo.cpp \
	geometry_lx.cpp \
	geometry_mix.cpp \
	geometry_Wsklx.cpp

#headers
geometry_includedir=$(includedir)/geometry
geometry_include_HEADERS= \
	geometry_eo.hpp \
	geometry_lx.hpp \
	geometry_mix.hpp \
	geometry_Wsklx.hpp

#compile or not virtual parallelization
if USE_VNODES
libgeometry_la_SOURCES+=geometry_vir.cpp
geometry_include_HEADERS+=geometry_vir.hpp
endif

