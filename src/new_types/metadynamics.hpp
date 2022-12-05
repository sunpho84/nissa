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
    
    int def_after() const
    {
      return 30;
    }
    
    int each;
    
    int def_each() const
    {
      return 1;
    }
    
    double coeff;
    
    double def_coeff() const
    {
      return 1.0;
    }
    
    double width;
    
    double def_width() const
    {
      return 1.0;
    }
    
    double barr;
    
    double def_barr() const
    {
      return 10.0;
    }
    
    double force_out;

    double def_force_out() const
    {
      return 100.0;
    }
    
    double well_tempering;
    
    double def_well_tempering() const
    {
      return 0.0;
    }
    
    double bend;
    
    double def_bend() const
    {
      return 0.0;
    }
    
    int def_ngrid() const
    {
      return 0;
    }
    
    int ngrid;
    storable_vector_t<double> grid;
    
    void update(const int& isweep,
		const double& Q);
    
    double compute_pot_der(const double& x) const;
    
    double compute_pot(const double& x) const;
    
    void save(const char *path) const;
    
    void load(const char *path);
    
    void draw_force(const char *force_path) const;
    
    void init();
    
    void read_pars();
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const;
    
    int is_nonstandard() const
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
