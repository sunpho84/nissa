#ifndef _METADYNAMICS_HPP
#define _METADYNAMICS_HPP

#include <math.h>
#include <string>
#include <vector>

#include "base/debug.hpp"
#include "io/ILDG_File.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  struct meta_pars_t
  {
    int after;
    int each;
    double coeff;
    double width;
    double barr;
    double force_out;
    double well_tempering;
    double bend;
    
    int def_after() {return 30;}
    int def_each() {return 1;}
    double def_coeff() {return 1.0;}
    double def_width() {return 1.0;}
    double def_barr() {return 10.0;}
    double def_force_out() {return 100.0;}
    double def_well_tempering() {return 0.0;}
    double def_bend() {return 0.0;}
    int def_ngrid() {return 0;}
    
    int ngrid;
    storable_vector_t<double> grid;
    
    void update(int isweep,double Q);
    double compute_pot_der(double x);
    double compute_pot(double x);
    void save(const char *path);
    void load(const char *path);
    void draw_force(const char *force_path);
    void init();
    void read_pars();
    
    int master_fprintf(FILE *fout,int full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false);
    
    int is_nonstandard()
    {
      return
	after!=def_after() or
	each!=def_each() or
	coeff!=def_coeff() or
	width!=def_width() or
	barr!=def_barr() or
	force_out!=def_force_out() or
	well_tempering!=def_well_tempering() or
	bend!=def_bend() or
	ngrid!=def_ngrid();
    }
    
    meta_pars_t() :
      after(def_after()),
      each(def_each()),
      coeff(def_coeff()),
      width(def_width()),
      barr(def_barr()),
      force_out(def_force_out()),
      well_tempering(def_well_tempering()),
      bend(def_bend()),
      ngrid(def_ngrid()){}
  };
}

#endif
