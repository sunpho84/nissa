#ifndef _ALL_RECTANGLES_HPP
#define _ALL_RECTANGLES_HPP

#include "operations/smearing/smooth.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //parameters to measure all rectangles path
  struct all_rects_meas_pars_t
  {
    int each;
    
    int after;
    
    std::string path;
    
    smooth_pars_t temp_smear_pars;
    
    smooth_pars_t spat_smear_pars;
    
    int Tmin,Tmax,Dmin,Dmax;
    
    int def_each() const
    {
      return 1;
    }
    
    int def_after() const
    {
      return 0;
    }
    
    std::string def_path() const
    {
      return "rectangles";
    }
    
    int def_Tmin() const
    {
      return 1;
    }
    
    int def_Tmax() const
    {
      return 9;
    }
    
    int def_Dmin() const
    {
      return 1;
    }
    
    int def_Dmax() const
    {
      return 9;
    }
    
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"MeasAllRects\n";
      if(each!=def_each() or full) os<<" Each\t\t=\t"<<each<<"\n";
      if(after!=def_after() or full) os<<" After\t\t=\t"<<after<<"\n";
      if(path!=def_path() or full) os<<" Path\t\t=\t\""<<path.c_str()<<"\"\n";
      if(Dmin!=def_Dmin() or full) os<<" Dmin\t\t=\t"<<Dmin<<"\n";
      if(Dmax!=def_Dmax() or full) os<<" Dmax\t\t=\t"<<Dmax<<"\n";
      if(Tmin!=def_Tmin() or full) os<<" Tmin\t\t=\t"<<Tmin<<"\n";
      if(Tmax!=def_Tmax() or full) os<<" Tmax\t\t=\t"<<Tmax<<"\n";
      if(spat_smear_pars.is_nonstandard() or full) os<<" Spatial "<<spat_smear_pars.get_str(full);
      if(temp_smear_pars.is_nonstandard() or full) os<<" Temporal "<<temp_smear_pars.get_str(full);
      
      return os.str();
    }
    
    int is_nonstandard() const
    {
      return
	each!=def_each() or
	after!=def_after() or
	path!=def_path() or
	Tmin!=def_Tmin() or
	Tmax!=def_Tmax() or
	Dmin!=def_Dmin() or
	Dmax!=def_Dmax() or
	temp_smear_pars.is_nonstandard() or
	spat_smear_pars.is_nonstandard();
    }
    
    all_rects_meas_pars_t() :
      each(def_each()),
      after(def_after()),
      path(def_path()),
      Tmin(def_Tmin()),
      Tmax(def_Tmax()),
      Dmin(def_Dmin()),
      Dmax(def_Dmax())
    {
    }
  };
  
  void measure_all_rectangular_paths(const all_rects_meas_pars_t& pars,
				     const LxField<quad_su3>& ori_conf,
				     const int& iconf,
				     const int& create_output_file);
  void measure_all_rectangular_paths_old(all_rects_meas_pars_t *pars,quad_su3  *conf,int iconf,int create_output_file);
  void measure_all_rectangular_paths(all_rects_meas_pars_t *pars,eo_ptr<quad_su3> conf,int iconf,int create_output_file);
}

#endif
