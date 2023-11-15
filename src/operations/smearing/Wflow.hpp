#ifndef _WFLOW_HPP
#define _WFLOW_HPP

#include <sstream>

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/backfield.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "operations/covariant_derivative.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  namespace Wflow
  {
    //we add with the new weight the previous one multiplied by the old weight
    inline void update_arg(LxField<quad_su3>& arg,
			   const LxField<quad_su3>& conf,
			   const double& dt,
			   const which_dir_t& dirs,
			   const int& iter)
    {
      conf.updateEdges();
      
      //Runge-Kutta coefficients
      constexpr std::array<double,3> RK_wn={1.0/4,8.0/9, 3.0/4};
      constexpr std::array<double,3> RK_wo={0,    -17.0/9, -1};
      
      //add the new argument of the exponential to the old one
      PAR(0,locVol,
	  CAPTURE(dirs,RK_wn,RK_wo,iter,dt,
		  TO_WRITE(arg),
		  TO_READ(conf)),
	  ivol,
	  {
	    for(int mu=0;mu<NDIM;mu++)
	      if(dirs[mu])
		{
		  //compute the new contribution
		  su3 staple,temp;
		  su3_put_to_zero(staple);
		  for(int inu=0;inu<NDIM-1;inu++)
		    {
		      int nu=perp_dir[mu][inu];
		      int A=ivol,B=loclxNeighup[A][nu],D=loclxNeighdw[A][nu],E=loclxNeighup[D][mu],F=loclxNeighup[A][mu];
		      unsafe_su3_prod_su3(       temp, conf[A][nu],conf[B][mu]);
		      su3_summ_the_prod_su3_dag(staple,temp,       conf[F][nu]);
		      unsafe_su3_dag_prod_su3(temp,    conf[D][nu],conf[D][mu]);
		      su3_summ_the_prod_su3(staple,    temp,       conf[E][nu]);
		    }
		  
		  //build Omega
		  su3 omega;
		  unsafe_su3_prod_su3_dag(omega,staple,conf[ivol][mu]);
		  
		  //compute Q and weight (the minus is there due to original stout)
		  su3 iQ,Q;
		  unsafe_su3_traceless_anti_hermitian_part(iQ,omega);
		  su3_prod_idouble(Q,iQ,-RK_wn[iter]*dt); //putting here the integration time
		  
		  //combine old and new
		  su3_prod_double(arg[ivol][mu],arg[ivol][mu],RK_wo[iter]);
		  su3_summassign(arg[ivol][mu],Q);
		}
	  });
    }
    
    //update the conf according to exp(i arg) conf_
    inline void update_conf(const LxField<quad_su3>& arg,
			    LxField<quad_su3>& conf,
			    const which_dir_t& dirs)
    {
      //integrate
      PAR(0,locVol,
	  CAPTURE(dirs,
		  TO_WRITE(conf),
		  TO_READ(arg)),
	  ivol,
	  {
	    for(int mu=0;mu<NDIM;mu++)
	      if(dirs[mu])
		{
		  su3 expiQ;
		  safe_hermitian_exact_i_exponentiate(expiQ,arg[ivol][mu]);
		  safe_su3_prod_su3(conf[ivol][mu],expiQ,conf[ivol][mu]);
		}
	  });
    }
    
    //call the 2-links Laplace operator, for staggered fields
    inline void Laplace_operator_switch(LxField<color0>&out,
					const LxField<quad_su3>& c,
					const which_dir_t& dirs,
					const LxField<color0>& in)
    {
       Laplace_operator_2_links(out,c,dirs,in);
    }
    
    //call the 1-link Laplace operator, for non-staggered fields
    inline void Laplace_operator_switch(LxField<spincolor>& out,
					const LxField<quad_su3>& c,
					const which_dir_t& dirs,
					const LxField<spincolor>& in)
    {
      Laplace_operator(out,c,dirs,in);
    }
  }
  
  //structure to Wilson flow
  struct Wflow_pars_t
  {
    int nflows;
    
    double dt;
    
    int nrecu;
    
    double def_nflows() const
    {
      return 50;
    }
    
    double def_dt() const
    {
      return 0.2;
    }
    
    int def_nrecu() const
    {
      return 5;
    }
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"WFlow\n";
      if(full or is_nonstandard())
	{
	  if(full or nflows!=def_nflows()) os<<" NFlows\t=\t"<<nflows<<"\n";
	  if(full or dt!=def_dt()) os<<" FlowStep\t=\t"<<dt<<"\n";
	  if(full or nrecu!=def_nrecu()) os<<" NRecu\t=\t"<<nrecu<<"\n";
	}
      
      return os.str();
    }
    
    int is_nonstandard() const
    {
      return
	nflows!=def_nflows() or
	dt!=def_dt() or
	nrecu!=def_nrecu();
    }
    
    Wflow_pars_t() :
      nflows(def_nflows()),
      dt(def_dt()),
      nrecu(def_nrecu())
    {
    }
  };
  
  /////////////////////////////////////////////////// fermions /////////////////////////////////////////////////////////
  
  //store conf used to perform fermionic flow
  template <typename T=color0,
	    int nint_steps=3>
  struct internal_fermion_flower_t
  {
    int nd;
    
    //time step
    double dt;
    
    //the steps
    std::vector<LxField<quad_su3>> conf;
    
    //dirs to smear
    which_dir_t dirs;
    
    //storage for staples
    LxField<quad_su3> arg;
    
    //creator
    internal_fermion_flower_t(const double& dt,
			      const which_dir_t& dirs) :
      nd(locVol*sizeof(T)/sizeof(double)),
      dt(dt),
      conf(nint_steps,"conf"),
      dirs(dirs),
      arg("arg")
    {
    }
    
    //add or remove backfield
    void add_or_rem_backfield_to_confs(const bool& add_rem,
				       const OldEoField<quad_u1>& u1)
    {
      for(int i=0;i<nint_steps;i++)
	add_or_rem_backfield_with_or_without_stagphases_to_conf(conf[i],add_rem,u1,true);
    }
    
    //setup from a given conf - nb: the conf is always evolved forward
    void generate_intermediate_steps(const LxField<quad_su3>& ori_conf)
    {
      //store the original conf in conf[0]
      conf[0]=ori_conf;
      
      //first two steps of the gluon R.K
      for(int iter=0;iter<nint_steps-1;iter++)
	{
	  conf[iter+1]=conf[iter];
	  Wflow::update_arg(arg,conf[iter+1],dt,dirs,iter);
	  Wflow::update_conf(arg,conf[iter+1],dirs);
	}
      
      // for(int iter=0;iter<nint_steps;iter++)
      // 	master_printf("plaquette internal steo %d: %.16lg\n",iter,global_plaquette_lx_conf(conf[iter]));
    }
    
    //destroyer
    ~internal_fermion_flower_t()
    {
      conf.clear();
    }
  };
  
  //forward fermion flower
  template <typename T,
	    int nint_steps=3>
  struct fermion_flower_t :
    public internal_fermion_flower_t<T,nint_steps>
  {
    //aux fields
    LxField<T> df0,df1,df2,f1,f2;
    
    //creator
    fermion_flower_t(double dt,
		     const which_dir_t& ext_dirs) :
      internal_fermion_flower_t<T,nint_steps>(dt,ext_dirs),
      df0("df0",WITH_HALO),
      df1("df1",WITH_HALO),
      df2("df2",WITH_HALO),
      f1("f1",WITH_HALO),
      f2("f2",WITH_HALO)
    {
    }
    
    //flow a field
    void flow_fermion(const LxField<T>& field)
    {
      auto& conf=this->conf;
      const double& dt=this->dt;
      const int& nd=this->nd;
      
      const LxField<T>& f0=field;
      LxField<T>& df0=this->df0;
      LxField<T>& f1=this->f1;
      
      //zero step: phi1 = phi0 + de0/4
      Wflow::Laplace_operator_switch(df0,conf[0],this->dirs,f0);
      
      FOR_EACH_SITE_DEG_OF_FIELD(f1,
				 CAPTURE(dt,
					 TO_WRITE(f1),
					 TO_READ(f0),
					 TO_READ(df0)),
				 site,deg,
				 {
				   f1(site,deg)=f0(site,deg)+df0(site,deg)*dt/4;
				 });
      
      crash("reimplement"); (void)nd;
      
      // //first step: phi2 = phi0 + de1*8/9 - df0*2/9
      // Wflow::Laplace_operator_switch(df1,conf[1],this->dirs,f1);
      // double_vector_summ_double_vector_prod_double((double*)f2,(double*)f0,(double*)df1,8.0*dt/9,nd);
      // double_vector_summassign_double_vector_prod_double((double*)f2,(double*)df0,-2.0*dt/9,nd);
      // //second step: f = phi3 = phi1 + df2*3/4
      // Wflow::Laplace_operator_switch(df2,conf[2],this->dirs,f2);
      // double_vector_summ_double_vector_prod_double((double*)f0,(double*)f1,(double*)df2,3.0*dt/4,nd);
    }
    
    //make the tail of the flow the head for next step
    void prepare_for_next_flow(LxField<quad_su3>& ext_conf)
    {
      if(nint_steps!=4) crash("not flown to last step!");
      ext_conf=this->conf[3];
    }
  };
  
  //adjoint fermion flower
  template <typename T=color0,
	    int nint_steps=3>
  struct fermion_adjoint_flower_t :
    public internal_fermion_flower_t<T,nint_steps>
  {
    //aux fields
    LxField<T> l1,l2;
    
    //creator
    fermion_adjoint_flower_t(double dt,const which_dir_t& ext_dirs) :
      internal_fermion_flower_t<T,nint_steps>(dt,ext_dirs),
      l1("l1",WITH_HALO),
      l2("l2",WITH_HALO)
    {
    }
    
    //flow a field
    void flow_fermion(LxField<T>& field)
    {
      crash("reimplement");
      
      // LxField<quad_su3> *conf=this->conf;
      // const double &dt=this->dt;
      // const int &nd=this->nd;
      
      // T *l3=field,*l0=l3;
      
      // //zero step: l2 = d2l3*3/4
      // Wflow::Laplace_operator_switch(l2,conf[2],this->dirs,l3);
      // double_vector_prodassign_double((double*)l2,3.0*dt/4,nd);
      // //first step: l1 = l3 + d1l2*8/9
      // Wflow::Laplace_operator_switch(l1,conf[1],this->dirs,l2);
      // double_vector_summ_double_vector_prod_double((double*)l1,(double*)l3,(double*)l1,8.0*dt/9,nd);
      // //second step: l0 = l1 + l2 + d0 (l1 - l2*8/9)/4
      // double_vector_summ((double*)l0,(double*)l1,(double*)l2,nd);                            //l0 = l1 + l2
      // double_vector_summassign_double_vector_prod_double((double*)l1,(double*)l2,-8.0/9,nd); //l1 = l1 - l2*8/9
      // Wflow::Laplace_operator_switch(l2,conf[0],this->dirs,l1);                              //l2 = d0 (l1 - l2*8/9)
      // double_vector_summassign_double_vector_prod_double((double*)l0,(double*)l2,dt/4,nd);   //l0+= d0 (l1 - l2*8/9)/4
    }
  };
  

  //flow for the given time for a dt using 1006.4518 appendix C
  inline void Wflow_lx_conf(LxField<quad_su3>& conf,
			    const double& dt,
			    const which_dir_t& dirs)
  {
    //storage for staples
    LxField<quad_su3> arg("arg");
    arg.reset();
    
    //we write the 4 terms of the Runge Kutta scheme iteratively
    for(int iter=0;iter<3;iter++)
      {
	Wflow::update_arg(arg,conf,dt,dirs,iter);
	Wflow::update_conf(arg,conf,dirs);
      }
  }
}

#endif
