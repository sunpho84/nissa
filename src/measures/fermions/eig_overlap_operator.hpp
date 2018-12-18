#ifndef _EIG_OVERLAP_OPERATOR_HPP
#define _EIG_OVERLAP_OPERATOR_HPP

/////////////////////////////  DOBBIAMO METTERE GLI INCLUDE ////////////////////


namespace nissa
{
  struct eig_overlap_meas_pars_t : base_fermionic_meas_t
  {
    int neigs;
    double eig_precision;
    std::string def_path(){return "uccello_piviere";}
    int def_neigs(){return 5;}
    double def_eig_precision(){return 1e-5;}
    double def_M(){return 1;}
    double def_minerr(){return 1e-7;}    

    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout, "%s", get_str().cstr());}
    std::string get_str(bool full=false);


    int is_nonstandard()
    {
      return
        base_fermionic_meas_t::is_nonstandard()||
        path!=def_path() or
        neigs!=def_neigs() or
        eig_precsion!=def_eig_precision() or
        M!=def_M() or
        minerr!=def_minerr();
    }
   	eig_overlap_meas_pars_t() :
    	base_fermionic_meas_t(), neigs(def_neigs()), eig_precision(def_eig_precision())
    	{path=def_path();}
    	virtual ~eig_overlap_meas_pars_t(){}
  };
void measure_eig_overlap(complex**,eigvec, quad_su3**,conf,complex*, D_ov_eig_val, double M, double,minerr, int,neigs, double, eig_oprecision)
}
#endif
