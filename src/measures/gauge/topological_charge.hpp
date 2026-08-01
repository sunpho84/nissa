#ifndef _TOPOLOGICAL_CHARGE_HPP
#define _TOPOLOGICAL_CHARGE_HPP

#include "operations/smearing/smooth.hpp"

namespace nissa
{
  //parameters to measure topology properties
  struct TopMeasPars
  {
    int each;
    
    int after;
    
    int measCorr;
    
    std::string path;
    
    std::string corrPath;
    
    bool timeMomRep;
    
    smooth_pars_t smoothPars;
    
    int def_each() const
    {
      return 1;
    }
    
    int def_after() const
    {
      return 0;
    }
    
    int def_meas_corr() const
    {
      return 0;
    }
    
    std::string def_path() const
    {
      return "Topo";
    }
    
    std::string def_corr_path() const
    {
      return "TopoCorr";
    }
    
    int def_time_mom_rep() const
    {
      return 0;
    }
    
    bool is_nonstandard() const
    {
      return
	each!=def_each() or
	after!=def_after() or
	path!=def_path() or
	corrPath!=def_corr_path() or
	timeMomRep!=def_time_mom_rep() or
	smoothPars.is_nonstandard();
    }
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"MeasTop\n";
      if(each!=def_each() or full) os<<" Each\t\t=\t"<<each<<"\n";
      if(after!=def_after() or full) os<<" After\t\t=\t"<<after<<"\n";
      if(path!=def_path() or full) os<<" Path\t\t=\t\""<<path.c_str()<<"\"\n";
      if(measCorr!=def_meas_corr() or full) os<<" MeasCorr\t=\t"<<measCorr<<"\n";
      if(corrPath!=def_corr_path() or full) os<<" CorrPath\t=\t\""<<corrPath<<"\"\n";
      if(timeMomRep!=def_time_mom_rep() or full) os<<" TimeMomRep\t=\t"<<timeMomRep<<"\n";
      os<<smoothPars.get_str(full);
      
      return os.str();
    }
    
    TopMeasPars() :
      each(def_each()),
      after(def_after()),
      measCorr(def_meas_corr()),
      path(def_path()),
      corrPath(def_corr_path()),
      timeMomRep(def_time_mom_rep())
    {
    }
  };
  
  void measure_topology_eo_conf(const TopMeasPars &pars,
				const EoField<quad_su3>& unsmoothed_conf_eo,
				const int& iconf,
				const bool& conf_created);
  
  void measure_topology_lx_conf(const TopMeasPars &pars,
				const LxField<quad_su3>& unsmoothed_conf,
				const int& iconf,
				const bool& conf_created,
				const bool& preserve_unsmoothed);
  
  void local_topological_charge(double *charge,quad_su3 *conf);
  
  double total_topological_charge_eo_conf(const EoField<quad_su3>& eo_conf);
  
  void topological_staples(LxField<quad_su3>& staples,const LxField<quad_su3>& conf);
  
  double total_topological_charge_lx_conf(const LxField<quad_su3>& conf);
  
  /////////////////////////////////////////////////////////////////
  
  //This will calculate the six independent components of
  //              2*a^2*ig*P_{mu,nu}
  //for a single point. Note that P_{mu,nu} is still not anti-symmetric
  //please ensure to have communicated the edges outside!
  /*
    ^                   C--<-- B --<--Y
    |                   |  2  | |  1  |
    n                   |     | |     |
    u                   D-->--\X/-->--A
    |                   D--<--/X\--<--A
    -----mu---->        |  3  | |  4  |
    .                   |     | |     |
    .                   E-->-- F -->--G
    in order to have the anti-symmetric part, use
    the routine inside "clover_term"
  */
  template <typename L,
	    typename U>
  CUDA_HOST_AND_DEVICE
  void four_leaves_point(L&& leaves_summ,
			 const U& conf,
			 const int& X)
  {
    // if(conf.edgesAreValid) CRASH("communicate edges externally");
    
    for(int mu=0;mu<NDIM;mu++)
      {
	int A=loclxNeighup[X][mu];
	int D=loclxNeighdw[X][mu];
        
	for(int nu=mu+1;nu<NDIM;nu++)
	  {
	    int munu=edge_numb[mu][nu];
	    
	    int B=loclxNeighup[X][nu];
	    int F=loclxNeighdw[X][nu];
            
	    int C=loclxNeighup[D][nu];
	    int E=loclxNeighdw[D][nu];
            
	    int G=loclxNeighdw[A][nu];
            
	    su3 temp1,temp2;
	    
	    //Leaf 1
	    unsafe_su3_prod_su3(temp1,conf[X][mu],conf[A][nu]);           //    B--<--Y
	    unsafe_su3_prod_su3_dag(temp2,temp1,conf[B][mu]);             //    |  1  |
	    unsafe_su3_prod_su3_dag(leaves_summ[munu],temp2,conf[X][nu]); //    |     |
	    /*                                                 */         //    X-->--A
            
	    //Leaf 2
	    unsafe_su3_prod_su3_dag(temp1,conf[X][nu],conf[C][mu]);       //    C--<--B
	    unsafe_su3_prod_su3_dag(temp2,temp1,conf[D][nu]);             //    |  2  |
	    unsafe_su3_prod_su3(temp1,temp2,conf[D][mu]);                 //    |     |
	    su3_summ(leaves_summ[munu],leaves_summ[munu],temp1);          //    D-->--X
	    
	    //Leaf 3
	    unsafe_su3_dag_prod_su3_dag(temp1,conf[D][mu],conf[E][nu]);    //   D--<--X
	    unsafe_su3_prod_su3(temp2,temp1,conf[E][mu]);                  //   |  3  |
	    unsafe_su3_prod_su3(temp1,temp2,conf[F][nu]);                  //   |     |
	    su3_summ(leaves_summ[munu],leaves_summ[munu],temp1);           //   E-->--F
            
	    //Leaf 4
	    unsafe_su3_dag_prod_su3(temp1,conf[F][nu],conf[F][mu]);         //  X--<--A
	    unsafe_su3_prod_su3(temp2,temp1,conf[G][nu]);                   //  |  4  |
	    unsafe_su3_prod_su3_dag(temp1,temp2,conf[X][mu]);               //  |     |
	    su3_summ(leaves_summ[munu],leaves_summ[munu],temp1);            //  F-->--G
            
	    munu++;
	  }
      }
  }
}

#endif
