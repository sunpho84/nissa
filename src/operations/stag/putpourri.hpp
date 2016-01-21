#ifndef _PUTPOURRI_HPP
#define _PUTPOURRI_HPP

namespace nissa
{
  struct fermionic_putpourri_meas_pars_t
  {
    int each;
    int after;
    std::string path;
    double residue;
    int compute_susc;
    int ncopies;
    int nhits;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    std::string def_path(){return "lavanda";}
    double def_residue(){return 1e-12;}
    int def_compute_susc(){return 0;}
    int def_ncopies(){return 1;}
    int def_nhits(){return 1;}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    int is_nonstandard()
    {
      return
	each!=def_each()||
	after!=def_after()||
	path!=def_path()||
	residue!=def_residue()||
	compute_susc!=def_compute_susc()||
	ncopies!=def_ncopies()||
	nhits!=def_nhits();
    }
    
    fermionic_putpourri_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      path(def_path()),
      residue(def_residue()),
      compute_susc(def_compute_susc()),
      ncopies(def_ncopies()),
      nhits(def_nhits()) {}
  };
  
  void measure_fermionic_putpourri(quad_su3 **conf,theory_pars_t &theory_pars,fermionic_putpourri_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
