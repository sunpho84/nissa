#ifndef _RENDENS_HPP
#define _RENDENS_HPP

namespace nissa
{
  struct quark_rendens_meas_pars_t
  {
    int each;
    int after;
    int max_order;
    std::string path;
    double residue;
    int ncopies;
    int nhits;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    int def_max_order(){return 2;}
    std::string def_path(){return "rende";}
    double def_residue(){return 1e-12;}
    int def_ncopies(){return 1;}
    int def_nhits(){return 1;}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    int is_nonstandard()
    {
      return
	each!=def_each()||
	after!=def_after()||
	max_order!=def_max_order()||
	path!=def_path()||
	residue!=def_residue()||
	ncopies!=def_ncopies()||
	nhits!=def_nhits();
    }
    
    quark_rendens_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      max_order(def_max_order()),
      path(def_path()),
      residue(def_residue()),
      ncopies(def_ncopies()),
      nhits(def_nhits()) {}
  };
  
  void measure_quark_rendens(quad_su3 **conf,theory_pars_t &theory_pars,quark_rendens_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
