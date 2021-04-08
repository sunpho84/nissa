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
    void update_arg(quad_su3 *arg,quad_su3 *conf,double dt,bool *dirs,int iter);
    void update_conf(quad_su3 *arg,quad_su3 *conf,bool *dirs);
    
    //call the 2-links Laplace operator, for staggered fields
    inline void Laplace_operator_switch(color *out,quad_su3 *c,bool *dirs,color *in)
    {
      Laplace_operator_2_links(out,c,dirs,in);
    }
    
    //call the 1-link Laplace operator, for non-staggered fields
    inline void Laplace_operator_switch(spincolor *out,quad_su3 *c,bool *dirs,spincolor *in)
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
    double def_nflows(){return 50;}
    double def_dt(){return 0.2;}
    int def_nrecu(){return 5;}
    
    int master_fprintf(FILE *fout,bool full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(bool full=false)
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
    
    int is_nonstandard()
    {
      return
	nflows!=def_nflows() or
	dt!=def_dt() or
	nrecu!=def_nrecu();
    }
    
    Wflow_pars_t() :
      nflows(def_nflows()),
      dt(def_dt()),
      nrecu(def_nrecu()) {}
  };
  
  /////////////////////////////////////////////////// fermions /////////////////////////////////////////////////////////
  
  //store conf used to perform fermionic flow
  template <typename T=color,
	    int nint_steps=3>
  struct internal_fermion_flower_t
  {
    int nd;
    
    //time step
    double dt;
    //the steps
    quad_su3 *conf[nint_steps];
    //dirs to smear
    bool dirs[NDIM];
    //storage for staples
    quad_su3 *arg;
    //creator
    internal_fermion_flower_t(double dt,bool *ext_dirs) : dt(dt)
    {
      //copy dirs
      for(int mu=0;mu<NDIM;mu++) dirs[mu]=ext_dirs[mu];
      //allocate confs
      for(int iter=0;iter<nint_steps;iter++)
	conf[iter]=nissa_malloc("conf",(locVol+bord_vol+edge_vol).nastyConvert(),quad_su3);
      //alllocate staple
      arg=nissa_malloc("arg",locVol.nastyConvert(),quad_su3);
      
      nd=locVol.nastyConvert()*sizeof(T)/sizeof(double);
    }
    
    //add or remove backfield
    void add_or_rem_backfield_to_confs(bool add_rem,eo_ptr<quad_u1> u1)
    {
      for(int i=0;i<nint_steps;i++) add_or_rem_backfield_with_or_without_stagphases_to_conf(conf[i],add_rem,u1,true);
    }
    
    //setup from a given conf - nb: the conf is always evolved forward
    void generate_intermediate_steps(quad_su3 *ori_conf)
    {
      //store the original conf in conf[0]
      vector_copy(conf[0],ori_conf);
      
      //first two steps of the gluon R.K
      for(int iter=0;iter<nint_steps-1;iter++)
	{
	  vector_copy(conf[iter+1],conf[iter]);
	  Wflow::update_arg(arg,conf[iter+1],dt,dirs,iter);
	  Wflow::update_conf(arg,conf[iter+1],dirs);
	}
      
      // for(int iter=0;iter<nint_steps;iter++)
      // 	master_printf("plaquette internal steo %d: %.16lg\n",iter,global_plaquette_lx_conf(conf[iter]));
    }
    
    //destroyer
    ~internal_fermion_flower_t()
    {
      for(int i=0;i<nint_steps;i++)
	nissa_free(conf[i]);
      nissa_free(arg);
    }
  };
  
  //forward fermion flower
  template <typename T,
	    int nint_steps=3>
  struct fermion_flower_t : public internal_fermion_flower_t<T,nint_steps>
  {
    //aux fields
    T *df0,*df1,*df2,*f1,*f2;
    
    //creator
    fermion_flower_t(double dt,bool *ext_dirs) : internal_fermion_flower_t<T,nint_steps>(dt,ext_dirs)
    {
      auto size=(locVol+bord_vol).nastyConvert();
      df0=nissa_malloc("df0",size,T);
      df1=nissa_malloc("df1",size,T);
      df2=nissa_malloc("df2",size,T);
      f1=nissa_malloc("f1",size,T);
      f2=nissa_malloc("f2",size,T);
    }
    
    //flow a field
    void flow_fermion(T *field)
    {
      quad_su3** conf=this->conf;
      double &dt=this->dt;
      int &nd=this->nd;
      
      T *f0=field;
      
      //zero step: phi1 = phi0 + de0/4
      Wflow::Laplace_operator_switch(df0,conf[0],this->dirs,f0);
      double_vector_summ_double_vector_prod_double((double*)f1,(double*)f0,(double*)df0,dt/4,nd);
      //first step: phi2 = phi0 + de1*8/9 - df0*2/9
      Wflow::Laplace_operator_switch(df1,conf[1],this->dirs,f1);
      double_vector_summ_double_vector_prod_double((double*)f2,(double*)f0,(double*)df1,8.0*dt/9,nd);
      double_vector_summassign_double_vector_prod_double((double*)f2,(double*)df0,-2.0*dt/9,nd);
      //second step: f = phi3 = phi1 + df2*3/4
      Wflow::Laplace_operator_switch(df2,conf[2],this->dirs,f2);
      double_vector_summ_double_vector_prod_double((double*)f0,(double*)f1,(double*)df2,3.0*dt/4,nd);
    }
    
    //make the tail of the flow the head for next step
    void prepare_for_next_flow(quad_su3 *ext_conf)
    {
      if(nint_steps!=4) crash("not flown to last step!");
      vector_copy(ext_conf,this->conf[3]);
    }
    
    //destroyer
    ~fermion_flower_t()
    {
      nissa_free(df0);
      nissa_free(df1);
      nissa_free(df2);
      nissa_free(f1);
      nissa_free(f2);
    }
  };
  
  //adjoint fermion flower
  template <typename T=color,
	    int nint_steps=3>
  struct fermion_adjoint_flower_t : public internal_fermion_flower_t<T,nint_steps>
  {
    //aux fields
    T *l1,*l2;
    
    //creator
    fermion_adjoint_flower_t(double dt,bool *ext_dirs) : internal_fermion_flower_t<T,nint_steps>(dt,ext_dirs)
    {
      auto size=(locVol+bord_vol).nastyConvert();
      
      l2=nissa_malloc("l2",size,T);
      l1=nissa_malloc("l1",size,T);
    }
    
    //flow a field
    void flow_fermion(T *field)
    {
      quad_su3 **conf=this->conf;
      double &dt=this->dt;
      int &nd=this->nd;
      
      T *l3=field,*l0=l3;
      
      //zero step: l2 = d2l3*3/4
      Wflow::Laplace_operator_switch(l2,conf[2],this->dirs,l3);
      double_vector_prodassign_double((double*)l2,3.0*dt/4,nd);
      //first step: l1 = l3 + d1l2*8/9
      Wflow::Laplace_operator_switch(l1,conf[1],this->dirs,l2);
      double_vector_summ_double_vector_prod_double((double*)l1,(double*)l3,(double*)l1,8.0*dt/9,nd);
      //second step: l0 = l1 + l2 + d0 (l1 - l2*8/9)/4
      double_vector_summ((double*)l0,(double*)l1,(double*)l2,nd);                            //l0 = l1 + l2
      double_vector_summassign_double_vector_prod_double((double*)l1,(double*)l2,-8.0/9,nd); //l1 = l1 - l2*8/9
      Wflow::Laplace_operator_switch(l2,conf[0],this->dirs,l1);                              //l2 = d0 (l1 - l2*8/9)
      double_vector_summassign_double_vector_prod_double((double*)l0,(double*)l2,dt/4,nd);   //l0+= d0 (l1 - l2*8/9)/4
    }
    
    //destroyer
    ~fermion_adjoint_flower_t()
    {
      nissa_free(l1);
      nissa_free(l2);
    }
  };
  
  void Wflow_lx_conf(quad_su3 *conf,double dt,bool *dirs=all_dirs);
}

#endif
