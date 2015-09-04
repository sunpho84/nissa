#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_mix.hpp"
#include "new_types/complex.hpp"
#include "new_types/float_128.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "topological_charge.hpp"

namespace nissa
{
  //This will calculate 2*a^2*ig*P_{mu,nu} for a single point
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
    the routine "Pmunu_term"
  */
  void four_leaves_point(as2t_su3 leaves_summ,quad_su3 *conf,int X)
  {
    if(!check_edges_valid(conf[0])) crash("communicate edges externally");
    
    int munu=0;
    for(int mu=0;mu<4;mu++)
      {
	int A=loclx_neighup[X][mu];
	int D=loclx_neighdw[X][mu];
        
	for(int nu=mu+1;nu<4;nu++)
	  {
	    int B=loclx_neighup[X][nu];
	    int F=loclx_neighdw[X][nu];
            
	    int C=loclx_neighup[D][nu];
	    int E=loclx_neighdw[D][nu];
            
	    int G=loclx_neighdw[A][nu];
            
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

  THREADABLE_FUNCTION_2ARG(four_leaves, as2t_su3*,Pmunu, quad_su3*,conf)
  {
    GET_THREAD_ID();
    communicate_lx_quad_su3_edges(conf);
    
    NISSA_PARALLEL_LOOP(X,0,loc_vol) four_leaves_point(Pmunu[X],conf,X);
    set_borders_invalid(Pmunu);
  }
  THREADABLE_FUNCTION_END

  //take anti-symmetric part and divide by 4
  void four_leaves_anti_symmetrize_fourth(as2t_su3 out,as2t_su3 in)
  {
    for(int munu=0;munu<6;munu++)
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  {
	    out[munu][ic1][ic2][0]=(in[munu][ic1][ic2][0]-in[munu][ic2][ic1][0])/4;
	    out[munu][ic1][ic2][1]=(in[munu][ic1][ic2][1]+in[munu][ic2][ic1][1])/4;
	  }
  }
  
  //takes the anti-symmetric part of the four-leaves
  THREADABLE_FUNCTION_2ARG(Pmunu_term, as2t_su3*,Pmunu,quad_su3*,conf)
  {
    GET_THREAD_ID();
    communicate_lx_quad_su3_edges(conf);
    
    //calculate U-U^dagger
    NISSA_PARALLEL_LOOP(X,0,loc_vol)
      {
	//compute the four leaves and anti-symmetrize
	as2t_su3 leaves_summ;
	four_leaves_point(leaves_summ,conf,X);
	four_leaves_anti_symmetrize_fourth(Pmunu[X],leaves_summ);
      }
    
    set_borders_invalid(Pmunu);
  }
  THREADABLE_FUNCTION_END
  
  //build the clover therm from anti-symmetrized four-leaves, following this decomposition
  /*
    (0,-P2+P3),     (-P1-P4,-P0+P5), (0,0),          (0,0)
    (P1+P4,-P0+P5), (0,P2-P3),       (0,0),          (0,0)
    (0,0),          (0,0),           (0,P2+P3),      (P1-P4,P0+P5)
    (0,0),          (0,0),           (-P1+P4,P0+P5), (0,-P2-P3)
    
    +iA   -B+iC  0      0
    +B+iC -iA    0      0
    0      0     +iD    -E+iF
    0      0     +E+iF  -iD
    
    A=P3-P2, B=P4+P1, C=P5-P0
    D=P3+P2, E=P4-P1, F=P5+P0
    
    +G  +H^+ 0    0
    +H  -G   0    0
    0    0   +I  +J^+
    0    0   +J  -I 
    
    out[0]=G=iA, out[1]=H=B+iC
    out[2]=I=iD, out[3]=J=E+iF
    
    NB: indeed Pi is anti-hermitian
  */
  void build_chromo_therm_from_anti_symmetric_four_leaves(quad_su3 out,as2t_su3 in)
  {
    su3 A,B,C,D,E,F;

    su3_subt(A,in[3],in[2]);
    su3_summ(B,in[4],in[1]);
    su3_subt(C,in[5],in[0]);
    su3_summ(D,in[3],in[2]);
    su3_subt(E,in[4],in[1]);
    su3_summ(F,in[5],in[0]);

    su3_prod_idouble(out[0],A,1);
    
    su3_copy(out[1],B);
    su3_summ_the_prod_idouble(out[1],C,1);
    
    su3_prod_idouble(out[2],D,1);
    
    su3_copy(out[3],E);
    su3_summ_the_prod_idouble(out[3],F,1);
  }
  
  //takes the anti-simmetric part of the four-leaves (optimized)
  THREADABLE_FUNCTION_2ARG(opt_Pmunu_term, quad_su3*,C,quad_su3*,conf)
  {
    GET_THREAD_ID();
    communicate_lx_quad_su3_edges(conf);
    
    //calculate U-U^dagger and store only the independent components
    NISSA_PARALLEL_LOOP(X,0,loc_vol)
      {
	//compute the four leaves and anti-symmetrize
	as2t_su3 leaves_summ;
	four_leaves_point(leaves_summ,conf,X);
	four_leaves_anti_symmetrize_fourth(leaves_summ,leaves_summ);
	build_chromo_therm_from_anti_symmetric_four_leaves(C[X],leaves_summ);
      }
    
    set_borders_invalid(C);
  }
  THREADABLE_FUNCTION_END
  
  //apply the chromo operator to the passed spinor site by site (not optimized)
  void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,as2t_su3 Pmunu,spincolor in)
  {
    for(int d1=0;d1<4;d1++)
      {
	color_put_to_zero(out[d1]);
	for(int imunu=0;imunu<6;imunu++)
	  {
	    color temp_d1;
	    unsafe_su3_prod_color(temp_d1,Pmunu[imunu],in[smunu_pos[d1][imunu]]);
	    for(int c=0;c<3;c++) complex_summ_the_prod(out[d1][c],smunu_entr[d1][imunu],temp_d1[c]);
	  }
      }
  }  
  void unsafe_apply_point_chromo_operator_to_spincolor_128(spincolor_128 out,as2t_su3 Pmunu,spincolor_128 in)
  {
    for(int d1=0;d1<4;d1++)
      {
	color_128_put_to_zero(out[d1]);
	for(int imunu=0;imunu<6;imunu++)
	  {
	    color_128 temp_d1;
	    unsafe_su3_prod_color_128(temp_d1,Pmunu[imunu],in[smunu_pos[d1][imunu]]);
	    for(int c=0;c<3;c++) complex_summ_the_64_prod_128(out[d1][c],smunu_entr[d1][imunu],temp_d1[c]);
	  }
      }
  }
  
  //apply the chromo operator to the passed spinor site by site (optimized)
  void unsafe_apply_opt_point_chromo_operator_to_spincolor(spincolor out,opt_as2t_su3 C,spincolor in)
  {
    unsafe_su3_prod_color(out[0],C[0],in[0]);
    su3_dag_summ_the_prod_color(out[0],C[1],in[1]);
    unsafe_su3_prod_color(out[1],C[1],in[0]);
    su3_subt_the_prod_color(out[1],C[0],in[1]);
    
    unsafe_su3_prod_color(out[2],C[2],in[2]);
    su3_dag_summ_the_prod_color(out[2],C[3],in[3]);
    unsafe_su3_prod_color(out[3],C[3],in[2]);
    su3_subt_the_prod_color(out[3],C[2],in[3]);
  }
  void unsafe_apply_opt_point_chromo_operator_to_spincolor_128(spincolor_128 out,opt_as2t_su3 C,spincolor_128 in)
  {
    unsafe_su3_prod_color_128(out[0],C[0],in[0]);
    su3_dag_summ_the_prod_color_128(out[0],C[1],in[1]);
    unsafe_su3_prod_color_128(out[1],C[1],in[0]);
    su3_subt_the_prod_color_128(out[1],C[0],in[1]);
    
    unsafe_su3_prod_color_128(out[2],C[2],in[2]);
    su3_dag_summ_the_prod_color_128(out[2],C[3],in[3]);
    unsafe_su3_prod_color_128(out[3],C[3],in[2]);
    su3_subt_the_prod_color_128(out[3],C[2],in[3]);
  }
  
  //apply the chromo operator to the passed spinor to the whole volume
  THREADABLE_FUNCTION_3ARG(unsafe_apply_chromo_operator_to_spincolor, spincolor*,out, as2t_su3*,Pmunu, spincolor*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      unsafe_apply_point_chromo_operator_to_spincolor(out[ivol],Pmunu[ivol],in[ivol]);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //apply the chromo operator to the passed spinor to the whole volume
  THREADABLE_FUNCTION_3ARG(unsafe_apply_chromo_operator_to_spincolor_128, spincolor_128*,out, as2t_su3*,Pmunu, spincolor_128*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      unsafe_apply_point_chromo_operator_to_spincolor_128(out[ivol],Pmunu[ivol],in[ivol]);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //apply the chromo operator to the passed spinor to the whole volume (the optimized way)
  THREADABLE_FUNCTION_3ARG(unsafe_apply_opt_chromo_operator_to_spincolor, spincolor*,out, opt_as2t_su3*,Cl, spincolor*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      unsafe_apply_opt_point_chromo_operator_to_spincolor(out[ivol],Cl[ivol],in[ivol]);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  THREADABLE_FUNCTION_3ARG(unsafe_apply_opt_chromo_operator_to_spincolor_128, spincolor_128*,out, opt_as2t_su3*,Cl, spincolor_128*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      unsafe_apply_opt_point_chromo_operator_to_spincolor_128(out[ivol],Cl[ivol],in[ivol]);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END

  //apply the chromo operator to the passed colorspinspin
  //normalization as in ape next
  THREADABLE_FUNCTION_3ARG(unsafe_apply_chromo_operator_to_colorspinspin, colorspinspin*,out, as2t_su3*,Pmunu, colorspinspin*,in)
  {
    spincolor temp1,temp2;
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      //Loop over the four source dirac indexes
      for(int id_source=0;id_source<4;id_source++) //dirac index of source
	{
	  //Switch the color_spinspin into the spincolor.
	  get_spincolor_from_colorspinspin(temp1,in[ivol],id_source);
	  
	  unsafe_apply_point_chromo_operator_to_spincolor(temp2,Pmunu[ivol],temp1);
	  
	  //Switch back the spincolor into the colorspinspin
	  put_spincolor_into_colorspinspin(out[ivol],temp2,id_source);
	}
    
    //invalidate borders
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END

  //apply the chromo operator to the passed su3spinspin
  //normalization as in ape next
  THREADABLE_FUNCTION_3ARG(unsafe_apply_chromo_operator_to_su3spinspin, su3spinspin*,out, as2t_su3*,Pmunu, su3spinspin*,in)
  {
    spincolor temp1,temp2;
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      //Loop over the four source dirac indexes
      for(int id_source=0;id_source<4;id_source++) //dirac index of source
	for(int ic_source=0;ic_source<3;ic_source++) //color index of source
	  {
	    //Switch the su3spinspin into the spincolor.
	    get_spincolor_from_su3spinspin(temp1,in[ivol],id_source,ic_source);
	    
	    unsafe_apply_point_chromo_operator_to_spincolor(temp2,Pmunu[ivol],temp1);
	    
	    //Switch back the spincolor into the colorspinspin
	    put_spincolor_into_su3spinspin(out[ivol],temp2,id_source,ic_source);
	  }
    
    //invalidate borders
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //measure the topological charge site by site
  THREADABLE_FUNCTION_2ARG(local_topological_charge, double*,charge, quad_su3*,conf)
  {
    double norm_fact=1/(128*M_PI*M_PI);
    
    as2t_su3 *leaves=nissa_malloc("leaves",loc_vol,as2t_su3);
    
    vector_reset(charge);
    
    //compute the clover-shape paths
    four_leaves(leaves,conf);
    
    //list the three combinations of plans
    int plan_id[3][2]={{0,5},{1,4},{2,3}};
    int sign[3]={1,-1,1};
    
    //loop on the three different combinations of plans
    GET_THREAD_ID();
    for(int iperm=0;iperm<3;iperm++)
      {
	//take the index of the two plans
	int ip0=plan_id[iperm][0];
	int ip1=plan_id[iperm][1];
	
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  {
	    //products
	    su3 clock,aclock;
	    unsafe_su3_prod_su3_dag(clock,leaves[ivol][ip0],leaves[ivol][ip1]);
	    unsafe_su3_prod_su3(aclock,leaves[ivol][ip0],leaves[ivol][ip1]);
	    
	    //take the trace
	    complex tclock,taclock;
	    su3_trace(tclock,clock);
	    su3_trace(taclock,aclock);
	    
	    //takes the combination with appropriate sign
	    charge[ivol]+=sign[iperm]*(tclock[RE]-taclock[RE])*norm_fact;
	  }
      }
    
    set_borders_invalid(charge);
    
    nissa_free(leaves);
  }
  THREADABLE_FUNCTION_END

  //total topological charge
  THREADABLE_FUNCTION_2ARG(total_topological_charge_lx_conf, double*,tot_charge, quad_su3*,conf)
  {
    GET_THREAD_ID();
    double *charge=nissa_malloc("charge",loc_vol,double);
    
    //compute local charge
    local_topological_charge(charge,conf);
    
    //summ over local volume
#ifndef REPRODUCIBLE_RUN
    double temp=0;
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      temp+=charge[ivol];
    
    *tot_charge=glb_reduce_double(temp);
#else
    //perform thread summ
    float_128 loc_thread_res={0,0};
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      float_128_summassign_64(loc_thread_res,charge[ivol]);
    
    float_128 temp;
    glb_reduce_float_128(temp,loc_thread_res);
    (*tot_charge)=temp[0];
#endif
    
    nissa_free(charge);
  }
  THREADABLE_FUNCTION_END

  //wrapper for eos case
  THREADABLE_FUNCTION_2ARG(total_topological_charge_eo_conf, double*,tot_charge, quad_su3**,eo_conf)
  {
    //convert to lx
    quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol+edge_vol,quad_su3);
    paste_eo_parts_into_lx_conf(lx_conf,eo_conf);
    
    total_topological_charge_lx_conf(tot_charge,lx_conf);
    
    nissa_free(lx_conf);
  }
  THREADABLE_FUNCTION_END

  //measure the topological charge
  void measure_topology_lx_conf(top_meas_pars_t &pars,quad_su3 *uncooled_conf,int iconf,bool conf_created,bool preserve_uncooled)
  {
    FILE *file=open_file(pars.path,conf_created?"w":"a");
    
    //allocate a temorary conf to be cooled
    quad_su3 *cooled_conf;
    if(preserve_uncooled)
      {
	cooled_conf=nissa_malloc("cooled_conf",loc_vol+bord_vol+edge_vol,quad_su3);
	vector_copy(cooled_conf,uncooled_conf);
      }
    else cooled_conf=uncooled_conf;
    
    //print curent measure and cool
    for(int istep=0;istep<=(pars.cool_pars.nsteps/pars.cool_pars.meas_each)*pars.cool_pars.meas_each;istep++)
      {
	if(istep%pars.cool_pars.meas_each==0)
	  {
	    double tot_charge;
	    total_topological_charge_lx_conf(&tot_charge,cooled_conf);
	    master_fprintf(file,"%d %d %16.16lg\n",iconf,istep,tot_charge);
	    verbosity_lv2_master_printf("Topological charge after %d cooling steps: %16.16lg, "
					"plaquette: %16.16lg\n",istep,tot_charge,global_plaquette_lx_conf(cooled_conf));
	  }
	if(istep!=pars.cool_pars.nsteps) cool_lx_conf(cooled_conf,pars.cool_pars.gauge_action,pars.cool_pars.overrelax_flag,pars.cool_pars.overrelax_exp);
      }
    
    //discard cooled conf
    if(preserve_uncooled) nissa_free(cooled_conf);
    
    close_file(file);
  }
  void measure_topology_eo_conf(top_meas_pars_t &pars,quad_su3 **uncooled_conf_eo,int iconf,bool conf_created)
  {
    quad_su3 *uncooled_conf_lx=nissa_malloc("uncooled_conf_lx",loc_vol+bord_vol+edge_vol,quad_su3);
    paste_eo_parts_into_lx_conf(uncooled_conf_lx,uncooled_conf_eo);
    measure_topology_lx_conf(pars,uncooled_conf_lx,iconf,conf_created,false);
    nissa_free(uncooled_conf_lx);
  }

  //compute the topological staples site by site
  THREADABLE_FUNCTION_2ARG(topological_staples, quad_su3*,staples, quad_su3*,conf)
  {
    GET_THREAD_ID();
    as2t_su3 *leaves=nissa_malloc("leaves",loc_vol+bord_vol+edge_vol,as2t_su3);
    
    //compute the clover-shape paths
    four_leaves(leaves,conf);
    //takes the anti-symmetric part (apart from a factor 2), in an horrendous wat
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int imunu=0;imunu<6;imunu++)
	{
	  color *u=leaves[ivol][imunu];
	  for(int ic1=0;ic1<3;ic1++)
	    for(int ic2=ic1;ic2<3;ic2++)
	      { //do not look here please, it is better to put a carpet on this uglyness
		u[ic2][ic1][0]=-(u[ic1][ic2][0]=u[ic1][ic2][0]-u[ic2][ic1][0]);
		u[ic2][ic1][1]=+(u[ic1][ic2][1]=u[ic1][ic2][1]+u[ic2][ic1][1]);
	      }
	}
    THREAD_BARRIER();
    set_borders_invalid(leaves);
    communicate_lx_as2t_su3_edges(leaves);
    
    //list the plan and coefficients for each staples
    int plan_perp[4][3]={{ 5, 4, 3},{ 5, 2, 1},{ 4, 2, 0},{ 3, 1, 0}};
    int plan_sign[4][3]={{+1,-1,+1},{-1,+1,-1},{+1,-1,+1},{-1,+1,-1}};
    
    //loop on the three different combinations of plans
    vector_reset(staples);
    NISSA_PARALLEL_LOOP(A,0,loc_vol)
      for(int mu=0;mu<4;mu++) //link direction
	for(int inu=0;inu<3;inu++)                   //  E---F---C   
	  {                                          //  |   |   | mu
	    int nu=perp_dir[mu][inu];                //  D---A---B   
	    //this gives the other pair element      //        nu    
	    int iplan=plan_perp[mu][inu];
	    
	    //takes neighbours
	    int B=loclx_neighup[A][nu];
	    int C=loclx_neighup[B][mu];
	    int D=loclx_neighdw[A][nu];
	    int E=loclx_neighup[D][mu];
	    int F=loclx_neighup[A][mu];
	    
	    //compute ABC, BCF and the full SU(3) staple, ABCF
	    su3 ABC,BCF,ABCF;
	    unsafe_su3_prod_su3(ABC,conf[A][nu],conf[B][mu]);
	    unsafe_su3_prod_su3_dag(BCF,conf[B][mu],conf[F][nu]);
	    unsafe_su3_prod_su3_dag(ABCF,ABC,conf[F][nu]);
	    
	    //compute ADE, DEF and the full SU(3) staple, ADEF
	    su3 ADE,DEF,ADEF;
	    unsafe_su3_dag_prod_su3(ADE,conf[D][nu],conf[D][mu]);
	    unsafe_su3_prod_su3(DEF,conf[D][mu],conf[E][nu]);
	    unsafe_su3_prod_su3(ADEF,ADE,conf[E][nu]);

	    //local summ and temp
	    su3 loc_staples,temp;
	    su3_put_to_zero(loc_staples);
	    //insert the leave in the four possible forward positions
	    
	    unsafe_su3_prod_su3(loc_staples,leaves[A][iplan],ABCF);     //insertion on A
	    unsafe_su3_prod_su3(temp,conf[A][nu],leaves[B][iplan]);
	    su3_summ_the_prod_su3(loc_staples,temp,BCF);                //insertion on B
	    unsafe_su3_prod_su3(temp,ABC,leaves[C][iplan]);
	    su3_summ_the_prod_su3_dag(loc_staples,temp,conf[F][nu]);    //insertion on C
	    su3_summ_the_prod_su3(loc_staples,ABCF,leaves[F][iplan]);   //insertion on F
	    
	    //insert the leave in the four possible backward positions
	    su3_summ_the_dag_prod_su3(loc_staples,leaves[A][iplan],ADEF);    //insertion on A
	    unsafe_su3_dag_prod_su3_dag(temp,conf[D][nu],leaves[D][iplan]);
	    su3_summ_the_prod_su3(loc_staples,temp,DEF);                     //insertion on D
	    unsafe_su3_prod_su3_dag(temp,ADE,leaves[E][iplan]);
	    su3_summ_the_prod_su3(loc_staples,temp,conf[E][nu]);             //insertion on E
	    su3_summ_the_prod_su3_dag(loc_staples,ADEF,leaves[F][iplan]);    //insertion on F
	    
	    //summ or subtract, according to the coefficient
	    if(plan_sign[mu][inu]==+1) su3_summassign(staples[A][mu],loc_staples);
	    else                       su3_subtassign(staples[A][mu],loc_staples);
	  }
    
    set_borders_invalid(staples);
    
    nissa_free(leaves);
  }
  THREADABLE_FUNCTION_END
  
  //store the topological charge if needed
  void topotential_pars_t::store_if_needed(quad_su3 **ext_conf,int iconf)
  {
    if(flag==2 && iconf%each==0 && iconf>=after)
      {
	double charge;
	quad_su3 *conf[2];
	if(stout_pars.nlev==0)
	  {
	    conf[0]=ext_conf[0];
	    conf[1]=ext_conf[1];
	  }
	else
	  {
	    conf[0]=nissa_malloc("stout_conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
	    conf[1]=nissa_malloc("stout_conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
	    stout_smear(conf,ext_conf,&stout_pars);
	  }
	
	//compute topocharge
	total_topological_charge_eo_conf(&charge,conf);
	update(iconf,charge);
	
	//free if needed
	if(stout_pars.nlev!=0)
	  {
	    nissa_free(conf[0]);
	    nissa_free(conf[1]);
	  }
      }
  }
}
