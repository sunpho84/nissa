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
    
    int master_fprintf(FILE *fout,int full) {return nissa::master_fprintf(fout,get_str().c_str());}
    std::string get_str(bool full=false);
    
    meta_pars_t() : after(30),each(1),coeff(1.0),width(1.0),barr(10.0),force_out(100.0),well_tempering(0.0),bend(0.0),ngrid(0){}
  };
}

#endif
