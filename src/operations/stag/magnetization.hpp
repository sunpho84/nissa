#ifndef _MAGNETIZATION_HPP
#define _MAGNETIZATION_HPP

namespace nissa
{
    struct magnetization_meas_pars_t
  {
    int each;
    int after;
    std::string path;
    double residue;
    int ncopies;
    int nhits;
    
    int def_each(){return 1;}
    int def_after(){return 0;}
    std::string def_path(){return "magnetization";}
    double def_residue(){return 1e-12;}
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
	ncopies!=def_ncopies()||
	nhits!=def_nhits();
    }
    
    magnetization_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      path(def_path()),
      residue(def_residue()),
      ncopies(def_ncopies()),
      nhits(def_nhits()) {}
  };
  
  void magnetization(complex *magn,quad_su3 **conf,quark_content_t *quark,color **rnd,color **chi,complex *point_magn,coords *arg,int mu,int nu);
  void magnetization(complex *magn,quad_su3 **conf,int quantization,quad_u1 **u1b,quark_content_t *quark,double residue,color **rnd);
  void measure_magnetization(quad_su3 **conf,theory_pars_t &theory_pars,magnetization_meas_pars_t &meas_pars,int iconf,int conf_created);
}

#endif
