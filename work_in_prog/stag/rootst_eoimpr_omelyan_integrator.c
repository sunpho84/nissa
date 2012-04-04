#pragma once

//unitarize the conf by explicitely inverting it
void eo_conf_unitarize_explicitely_inverting(quad_su3 **conf)
{
  addrem_stagphases_to_eo_conf(conf);
  
  for(int eo=0;eo<2;eo++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  su3_unitarize_explicitly_inverting(conf[eo][ivol][mu],conf[eo][ivol][mu]);
      
      set_borders_invalid(conf[eo]);
    }
      
  addrem_stagphases_to_eo_conf(conf);
}

// Evolve momenta according to the rooted staggered force
// calculate H=H-F*dt to evolve link momenta
// i.e calculate v(t+dt)=v(t)+a*dt
void evolve_momenta_with_rootst_force(quad_su3 **H,quad_su3 **conf,double beta,int nfl,quad_u1 ***u1b,color **pf,rat_approx *appr,double residue,double dt)
{
  //allocate force
  quad_su3 *F[2]={nissa_malloc("F0",loc_volh,quad_su3),nissa_malloc("F1",loc_volh,quad_su3)};
  
  //compute the force
  full_rootst_eoimpr_force(F,conf,beta,nfl,u1b,pf,appr,residue);
  
  //evolve
  master_printf("Evolve momenta with force, dt=%lg\n",dt);
  
  for(int eo=0;eo<2;eo++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  for(int ic1=0;ic1<3;ic1++)
	    for(int ic2=0;ic2<3;ic2++)
	      complex_subt_the_prod_idouble(H[eo][ivol][mu][ic1][ic2],F[eo][ivol][mu][ic1][ic2],dt);
      
      nissa_free(F[eo]);
    }
}

//eolve the configuration by using the computed momenta
void evolve_conf_with_momenta(quad_su3 **eo_conf,quad_su3 **H,double dt)
{
  //communitcate momenta borders
  master_printf("Evolving conf with momenta\n");
  
  //evolve
  for(int eo=0;eo<2;eo++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  {
	    su3 t1,t2;
	    su3_prod_with_idouble(t1,H[eo][ivol][mu],dt);
	    unsafe_su3_taylor_exponentiate(t2,t1,6);
	    safe_su3_prod_su3(eo_conf[eo][ivol][mu],t2,eo_conf[eo][ivol][mu]);
	  }
      set_borders_invalid(eo_conf[eo]);
    }
  
  master_printf("plaquette after: %.18lg\n",-global_plaquette_eo_conf(eo_conf));
}

// Omelyan integrator(cond-mat/0110438v1) for rooted staggered theory
// in summary, at each step it should make:
//
//     v1 = v(t) + a[r(t)]*lambda*dt
//     r1 = r(t) + v1*dt/2
//     v2 = v1 + a[r1]*(1 -2*lambda)*dt
//     r(t + dt) = r1 + v2*dt/2
//     v(t + h) = v2 + a[r(t + dt)]*lambda*dt.
//
// But it is sligthly optimized putting together first and last operations.
// This requires to distinguish first and last step from others:
//
//     only 1st step:
//      v1 = v(t) + a[r(t)]*lambda*dt
//     all the steps:
//      r1 = r(t) + v1*dt/2
//      v2 = v1 + a[r1]*(1 -2*lambda)*dt
//      r(t + dt) = r1 + v2*dt/2
//    not last step:
//     v1 = v2 + a[r(t + dt)]*2*lambda*dt
//    only last step:
//     v(t + h) = v2 + a[r(t + dt)]*lambda*dt
void omelyan_rootst_eoimpr_evolver(quad_su3 **H,quad_su3 **conf,double beta,int nfl,quad_u1 ***u1b,color **pf,rat_approx *appr,double residue,double traj_length,int nstep)
{
  //define step length ant its multiples
  const double lambda=0.1931833;
  double dt=traj_length/nstep,dth=dt/2;
  double ldt=dt*lambda,l2dt=2*lambda*dt,uml2dt=(1-2*lambda)*dt;  
  
  //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
  evolve_momenta_with_rootst_force(H,conf,beta,nfl,u1b,pf,appr,residue,ldt);
  
  //         Main loop
  for(int istep=0;istep<nstep;istep++)
    {
      master_printf("Starting step %d/%d\n",istep+1,nstep);
      
      //     Compute U(t+dt/2) i.e. r1=r(t)+v1*dt/2
      evolve_conf_with_momenta(conf,H,dth);
      //     Compute H(t+(1-2*lambda)*dt) i.e. v2=v1+a[r1]*(1-2*lambda)*dt
      evolve_momenta_with_rootst_force(H,conf,beta,nfl,u1b,pf,appr,residue,uml2dt);
      //     Compute U(t+dt/2) i.e. r(t+dt)=r1+v2*dt/2
      evolve_conf_with_momenta(conf,H,dth);
      //     Compute H(t+dt) i.e. v(t+dt)=v2+a[r(t+dt)]*lambda*dt !only! at last step
      if(istep==(nstep-1)) evolve_momenta_with_rootst_force(H,conf,beta,nfl,u1b,pf,appr,residue,ldt);
      //     Compute H(t+dt) i.e. v(t+dt)=v2+a[r(t+dt)]*2*lambda*dt at all but last step
      else                  evolve_momenta_with_rootst_force(H,conf,beta,nfl,u1b,pf,appr,residue,l2dt);
      
      //normalize the configuration
      eo_conf_unitarize_explicitely_inverting(conf);
    }
}
