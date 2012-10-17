#include <string.h>
#include <math.h>

#include "../new_types/new_types_definitions.h"
#include "../new_types/complex.h"
#include "../new_types/su3.h"
#include "../new_types/spin.h"
#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../base/routines.h"
#include "../base/communicate.h"
#include "../base/debug.h"
#include "../geometry/geometry_lx.h"
#include "../geometry/geometry_mix.h"

////////////////////////// Complicated things /////////////////////

//this is the general framework where to compute paths along the conf
//there are special routines are for optimized tasks

//init an su3 path
su3_path* su3_path_init(int ivol)
{
  su3_path *out=(su3_path*)malloc(sizeof(su3_path));
  out->ivol=ivol;
  su3_put_to_id(out->data);
  out->movements=0;
  out->next=NULL;
  
  return out;
}

//append to existing list of path
void su3_path_append(su3_path *out,su3_path *in)
{
  if(in->next!=NULL) crash("in->next!=NULL");
  if(out->next!=NULL) crash("out->next!=NULL");
  
  out->next=in;
}

//get the length of path
int su3_path_get_length(su3_path *in)
{return (in->movements)>>28;}

//get all the movements
int su3_path_get_movements(su3_path *in)
{return (in->movements)&268435455;}

//add a movement to su3_path
void su3_path_push_movement(su3_path *in,int mov)
{
  if(mov<0||mov>=8) crash("mov=%d, must be between 0 an 7",mov);
  
  int len=su3_path_get_length(in);
  int movements=mov<<(3*len)|su3_path_get_movements(in);
  len++;
  
  if(len>8) crash("too long path!");
  
  in->movements=(len<<28)|movements;
}

//add a movement to su3_path
void su3_path_pop_movement(su3_path *in)
{
  int len=su3_path_get_length(in)-1;
  if(len==-1) crash("empty list of movements!");
  
  int movements=su3_path_get_movements(in)>>3;
  
  in->movements=(len<<28)|movements;
}

//check if su3_path is terminated
int su3_path_check_finished(su3_path *in)
{return su3_path_get_length(in)==0;}

//take current movement
int su3_path_get_current_movement(su3_path *in)
{
  if(su3_path_check_finished(in)) crash("finished path!");
  
  return su3_path_get_movements(in)&7;
}

//square (the proto-plaquette)
/*
     
  C------D
n |      | 
u |      | 
  A--mu--B
  
The square path P_{mu,nu} is defined as U(A,mu)U(B,nu)U^(C,mu)U^(A,nu)=
=U(A,mu)U(B,nu)(U(A,nu)U^(C,mu))^=U(AB,munu)*U^(AC,numu)
*/

void squared_path(su3 square,quad_su3 *conf,int A,int mu,int nu)
{
  int B=loclx_neighup[A][mu];
  int C=loclx_neighup[A][nu];

  su3 ABD,ACD;

  unsafe_su3_prod_su3(ABD,conf[A][mu],conf[B][nu]);
  unsafe_su3_prod_su3(ACD,conf[A][nu],conf[C][mu]);
  unsafe_su3_prod_su3_dag(square,ABD,ACD);
}

//This calculate the global plaquette. It's not done in a very
//efficient way, but it's ok for our scope.
double global_plaquette_lx_conf(quad_su3 *conf)
{
  communicate_lx_quad_su3_borders(conf);
  
  su3 square;
  complex pl;
  double totplaq=0;
  nissa_loc_vol_loop(ivol)
    for(int idir=0;idir<4;idir++)
      for(int jdir=idir+1;jdir<4;jdir++)
	{
	  squared_path(square,conf,ivol,idir,jdir);
	  su3_trace(pl,square);
	  totplaq+=pl[0];
	}
  
  return glb_reduce_double(totplaq)/glb_vol/18;
}

//This calculate the variance of the global plaquette.
double global_plaquette_variance_lx_conf(quad_su3 *conf)
{
  communicate_lx_quad_su3_borders(conf);
  
  su3 square;
  complex pl;
  double totlocplaq=0,totlocplaq2=0;
  nissa_loc_vol_loop(ivol)
    for(int idir=0;idir<4;idir++)
      for(int jdir=idir+1;jdir<4;jdir++)
        {
          squared_path(square,conf,ivol,idir,jdir);
          su3_trace(pl,square);
          totlocplaq+=pl[0]/3;
          totlocplaq2+=(pl[0]/3)*(pl[0]/3);
        }
  
  double totplaq=glb_reduce_double(totlocplaq)/(6*glb_vol);
  double totplaq2=glb_reduce_double(totlocplaq2)/(6*glb_vol)-totplaq*totplaq;
  
  return sqrt(totplaq2);
}

/* compute the global plaquette on a e/o split conf, in a more efficient way
     
  C------D
n |      | 
u |      | 
  A--mu--B

  the temporal and spatial plaquette are computed separately
*/
void global_plaquette_eo_conf(double *totplaq,quad_su3 **conf)
{
  communicate_eo_quad_su3_borders(conf);
  
  double locplaq[2]={0,0};
  
  nissa_loc_volh_loop(A)
    for(int par=0;par<2;par++)
      for(int mu=0;mu<4;mu++)
	{
	  for(int nu=mu+1;nu<4;nu++)
	    {
	      //ACD and ABD path
	      su3 ABD,ACD;
	      unsafe_su3_prod_su3(ABD,conf[par][A][mu],conf[!par][loceo_neighup[par][A][mu]][nu]);
	      unsafe_su3_prod_su3(ACD,conf[par][A][nu],conf[!par][loceo_neighup[par][A][nu]][mu]);
	      
	      //compute tr(ABDC)
	      double tr=real_part_of_trace_su3_prod_su3_dag(ABD,ACD);
	      if(mu==0) locplaq[0]+=tr;
	      else      locplaq[1]+=tr;
	    }
	}
  
  //reduce double[2] as complex
  glb_reduce_complex(totplaq,locplaq);
  
  //normalize
  for(int ts=0;ts<2;ts++) totplaq[ts]/=glb_vol*3*3;
}

//return the average between spatial and temporary plaquette
double global_plaquette_eo_conf(quad_su3 **conf)
{
  double plaq[2];
  
  global_plaquette_eo_conf(plaq,conf);
  
  return (plaq[0]+plaq[1])/2;
}

//compute the staples along a particular dir, for a single site
void compute_point_staples_eo_conf_single_dir(su3 staple,quad_su3 **eo_conf,int A,int mu)
{
  communicate_eo_quad_su3_edges(eo_conf);
  
  su3_put_to_zero(staple);
  
  su3 temp1,temp2;
  for(int nu=0;nu<4;nu++)                   //  E---F---C   
    if(nu!=mu)                              //  |   |   | mu
      {                                     //  D---A---B   
	int p=loclx_parity[A];              //        nu    
	int B=loclx_neighup[A][nu];
	int F=loclx_neighup[A][mu];
	unsafe_su3_prod_su3(    temp1,eo_conf[p][loceo_of_loclx[A]][nu],eo_conf[!p][loceo_of_loclx[B]][mu]);
	unsafe_su3_prod_su3_dag(temp2,temp1,                            eo_conf[!p][loceo_of_loclx[F]][nu]);
	su3_summ(staple,staple,temp2);
	
	int D=loclx_neighdw[A][nu];
	int E=loclx_neighup[D][mu];
	unsafe_su3_dag_prod_su3(temp1,eo_conf[!p][loceo_of_loclx[D]][nu],eo_conf[!p][loceo_of_loclx[D]][mu]);
	unsafe_su3_prod_su3(    temp2,temp1,                             eo_conf[ p][loceo_of_loclx[E]][nu]);
	su3_summ(staple,staple,temp2);
      }
}

//compute the staples along all the four dirs
void compute_point_staples_eo_conf(quad_su3 staple,quad_su3 **eo_conf,int A)
{for(int mu=0;mu<4;mu++) compute_point_staples_eo_conf_single_dir(staple[mu],eo_conf,A,mu);}

//Compute rectangle staples. The routine loops over A and compute staples for 
//neighbouring points. If these are not on the local volume, they must be sent to the 
void compute_rectangle_staples_eo_conf(quad_su3 **staple,quad_su3 **eo_conf)
{
  communicate_eo_quad_su3_edges(eo_conf);
  
  //reset the staples - also borders are resetted
  for(int eo=0;eo<2;eo++) vector_reset(staple[eo]);
  
  nissa_loc_vol_loop(Alx)
    for(int mu=0;mu<4;mu++)
      for(int nu=mu+1;nu<4;nu++)
	{
	  int p=loclx_parity[Alx];
	  int A=loceo_of_loclx[Alx];
	  int B=loceo_neighup[p][A][nu];
	  int D=loceo_neighdw[p][A][nu];
	  int E=loceo_neighdw[!p][D][mu];
	  int F=loceo_neighup[p][A][mu];
	  
	  //compute DAB
	  su3 DAB;
	  unsafe_su3_prod_su3(DAB,eo_conf[!p][D][nu],eo_conf[p][A][nu]);
	  
	  //compute DABC
	  su3 DABC;
	  unsafe_su3_prod_su3(DABC,DAB,eo_conf[!p][B][mu]);
	  
	  //compute EDABC
	  su3 EDABC;
	  unsafe_su3_dag_prod_su3(EDABC,eo_conf[!p][D][mu],DABC);
	  
	  //compute EFC
	  su3 EFC;
	  unsafe_su3_prod_su3(EFC,eo_conf[p][E][nu],eo_conf[!p][F][nu]);
	  
	  //compute DEFC
	  su3 DEFC;
	  unsafe_su3_prod_su3(DEFC,eo_conf[!p][D][mu],EFC);
	  
	  //compute DEFCB
	  su3 DEFCB;
	  unsafe_su3_prod_su3_dag(DEFCB,DEFC,eo_conf[!p][B][mu]);
	  
	  // first of all the 2 staples when we consider the "horizontal" rectangle
	  //
	  //  E---F---C   
	  //  |   |   | mu
	  //  D---A---B   
	  //        nu    
	  
	  //closing with BC gives DABCFE staple     //  E-<-F-<-C
	  su3 DABCFE;				    //          |
	  unsafe_su3_prod_su3_dag(DABCFE,DABC,EFC); //  D->-A->-B
	  su3_summassign(staple[!p][D][mu],DABCFE);
	  
	  //closing with DE gives BADEFC staple     //  E->-F->-C
	  su3 BADEFC;				    //  |        
	  unsafe_su3_dag_prod_su3(BADEFC,DAB,DEFC); //  D-<-A-<-B
	  //su3_summassign(staple[!p][B][mu],BADEFC);
	  
	  // then the 4 staples when we consider the "vertical" one
	  //
	  //   B---C
	  //   |   |
	  //   A---F
	  // nu|   | 
	  //   D---E
	  //     mu
	  
	  //closing with DA gives ADEFCB staple
	  su3 ADEFCB;
	  unsafe_su3_dag_prod_su3(ADEFCB,eo_conf[!p][D][nu],DEFCB);
	  //su3_summassign(staple[p][A][nu],ADEFCB);
	    
	  //closing with AB gives DEFCBA staple
	  su3 DEFCBA;
	  unsafe_su3_prod_su3_dag(DEFCBA,DEFCB,eo_conf[p][A][nu]);
	  //su3_summassign(staple[!p][D][nu],DEFCBA);
	  
	  //closing with EF gives FEDABC staple
	  su3 FEDABC;
	  unsafe_su3_dag_prod_su3(FEDABC,eo_conf[p][E][nu],EDABC);
	  //su3_summassign(staple[!p][F][nu],FEDABC);
	  
	  //closing with FC gives EDABCF staple
	  su3 EDABCF;
	  unsafe_su3_prod_su3_dag(EDABCF,EDABC,eo_conf[!p][F][nu]);
	  //su3_summassign(staple[p][E][nu],EDABCF);
	}
}

//shift an su3 vector of a single step along the mu axis, in the positive or negative dir
void su3_vec_single_shift(su3 *u,int mu,int sign)
{
  //communicate borders
  communicate_lx_su3_borders(u);
  
  //choose the orthogonal dirs
  int nu=0,rho=0,sigma=0;
  switch(mu)
    {
    case 0:
      nu=1;
      rho=2;
      sigma=3;
      break;
    case 1:
      nu=0;
      rho=2;
      sigma=3;
      break;
    case 2:
      nu=0;
      rho=1;
      sigma=3;
      break;
    case 3:
      nu=0;
      rho=1;
      sigma=2;
      break;
    default:
      crash("mu>3 or mu<0");
    }
  
  //choose start, end and step
  int sh=(sign>0) ? -1 : +1;
  int st=(sign>0) ? loc_size[mu]-1 : 0;
  int en=(sign>0) ? 0 : loc_size[mu]-1 ;
  
  //loop over orthogonal dirs
  coords x;
  for(x[nu]=0;x[nu]<loc_size[nu];x[nu]++)
    for(x[rho]=0;x[rho]<loc_size[rho];x[rho]++)
      for(x[sigma]=0;x[sigma]<loc_size[sigma];x[sigma]++)
	{
	  //loop along the dir
	  x[mu]=st;
	  
	  //make a buffer in the case in which this dir is not parallelized
	  su3 buf;
	  int isite=loclx_of_coord(x);
	  if(nrank_dir[mu]==1)
	    su3_copy(buf,u[isite]);
	  
	  //loop on remaining dir
	  do
	    {
	      //find the site and the neighbour
	      int ineigh=(sign>0) ? loclx_neighdw[isite][mu] : loclx_neighup[isite][mu]; 
	      
	      //copy theneighbour in the site
	      su3_copy(u[isite],u[ineigh]);
	      	      
	      //advance
	      x[mu]+=sh;
	      if(x[mu]!=en+sh) isite=ineigh;
	    }
	  while(x[mu]!=en+sh);
	  
	  //if dir not parallelized, restore end
	  if(nrank_dir[mu]==1)
	    su3_copy(u[isite],buf);
	}
  
  //invalidate borders
  set_borders_invalid(u);
}

//compute the real part of the trace of the rectangle of size nstep_mu X nstep_nu in th mu|nu plane
void average_trace_of_rectangle_path(complex tra,quad_su3 *conf,int mu,int nu,int nstep_mu,int nstep_nu,su3 *u)
{
  communicate_lx_quad_su3_borders(conf);
  
  //reset the link product
  nissa_loc_vol_loop(ivol)
    su3_put_to_id(u[ivol]);

 //move along +mu
  for(int i=0;i<nstep_mu;i++)
    {
      //take the product
      nissa_loc_vol_loop(ivol)
	safe_su3_prod_su3(u[ivol],u[ivol],conf[ivol][mu]);
      
      su3_vec_single_shift(u,mu,+1);
    }

  //move along +nu
  for(int i=0;i<nstep_nu;i++)
    {
      //take the product
      nissa_loc_vol_loop(ivol)
	safe_su3_prod_su3(u[ivol],u[ivol],conf[ivol][nu]);
      
      su3_vec_single_shift(u,nu,+1);
    }
  
  //move along -mu
  for(int i=0;i<nstep_mu;i++)
    {
      su3_vec_single_shift(u,mu,-1);
      
      //take the product
      nissa_loc_vol_loop(ivol)
	safe_su3_prod_su3_dag(u[ivol],u[ivol],conf[ivol][mu]);
    }

  //move along -nu
  for(int i=0;i<nstep_nu;i++)
    {
      su3_vec_single_shift(u,nu,-1);
      
      //take the product
      nissa_loc_vol_loop(ivol)
	safe_su3_prod_su3_dag(u[ivol],u[ivol],conf[ivol][nu]);
    }
  
  //compute the trace
  complex loc_tra={0,0};
  nissa_loc_vol_loop(ivol) su3_summ_the_trace(loc_tra,u[ivol]);
  glb_reduce_complex(tra,loc_tra);
  
  //normalize
  complex_prod_double(tra,tra,1.0/glb_vol);
}

//return only the real part
double average_real_part_of_trace_of_rectangle_path(quad_su3 *conf,int mu,int nu,int nstep_mu,int nstep_nu,su3 *u)
{
  complex tra;
  average_trace_of_rectangle_path(tra,conf,mu,nu,nstep_mu,nstep_nu,u);
  
  return tra[0];
}

//compute the polyakov loop
void average_polyakov_loop(complex tra,quad_su3 *conf,int mu)
{
  su3 *u=nissa_malloc("u",loc_vol+bord_vol,su3);
  
  communicate_lx_quad_su3_borders(conf);
  
  //reset the link product
  nissa_loc_vol_loop(ivol)
    su3_put_to_id(u[ivol]);

  //move along +mu
  for(int i=0;i<glb_size[mu];i++)
    {
      //take the product
      nissa_loc_vol_loop(ivol)
	safe_su3_prod_su3(u[ivol],u[ivol],conf[ivol][mu]);
      
      su3_vec_single_shift(u,mu,+1);
    }
  
  //compute the trace; since we reduce over all the volume there are glb_size[mu] replica
  complex loc_tra={0,0};
  nissa_loc_vol_loop(ivol)
      su3_summ_the_trace(loc_tra,u[ivol]);
  glb_reduce_complex(tra,loc_tra);
  complex_prodassign_double(tra,1.0/glb_vol/3.0);
  
  nissa_free(u);
}

//definition in case of eos conf
void average_polyakov_loop_of_eos_conf(complex tra,quad_su3 **eo_conf,int mu)
{
  quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol,quad_su3);
  paste_eo_parts_into_lx_conf(lx_conf,eo_conf);
  
  average_polyakov_loop(tra,lx_conf,mu);
  
  nissa_free(lx_conf);
}

void Pline(su3 *Pline,quad_su3 *conf)
{
  communicate_lx_quad_su3_borders(conf);
  
  int X,iX[4],T=loc_size[0];
  
  su3 *U0=nissa_malloc("U0",loc_size[0]+1,su3);
  su3 *plf=nissa_malloc("plf",loc_size[0]+1,su3);
  su3 *plb=nissa_malloc("plb",loc_size[0]+1,su3);
  
  MPI_Status status_f,status_b;
  int tag_f=0,tag_b=1;
  su3 aux,plf_dw_buffer,plb_up_buffer;
  
  //local loop
  for(iX[1]=0;iX[1]<loc_size[1];iX[1]++)
    for(iX[2]=0;iX[2]<loc_size[2];iX[2]++)
      for(iX[3]=0;iX[3]<loc_size[3];iX[3]++)
	{
	  //forward line
	  for(iX[0]=0;iX[0]<T;iX[0]++)
	    {
	      X=loclx_of_coord(iX);
	      unsafe_su3_hermitian(U0[iX[0]+1],conf[X][0]);
             }
            su3_copy(plf[1],U0[1]);
            for(int t=2;t<=T;t++) unsafe_su3_prod_su3(plf[t],U0[t],plf[t-1]);

	    if(rank_coord[0]<nrank_dir[0]-1)
	      {
		MPI_Send(plf[T],1,MPI_SU3,rank_neighup[0],tag_f,MPI_COMM_WORLD);
		MPI_Recv(plb_up_buffer,1,MPI_SU3,rank_neighup[0],tag_b,MPI_COMM_WORLD,&status_b);
	      }

	    //backward line
            for(iX[0]=0;iX[0]<T;iX[0]++)
	      {
		X=loclx_of_coord(iX);
		su3_copy(U0[iX[0]+1],conf[X][0]);
	      }
            su3_copy(plb[T],U0[T]);
            for(int t=1;t<T;t++) unsafe_su3_prod_su3(plb[T-t],U0[T-t],plb[T-t+1]);

	   if(rank_coord[0]>0)
	     {
	       MPI_Send(plb[1],1,MPI_SU3,rank_neighdw[0],tag_b,MPI_COMM_WORLD);
	       MPI_Recv(plf_dw_buffer,1,MPI_SU3,rank_neighdw[0],tag_f,MPI_COMM_WORLD,&status_f);
	     }

	    //forward P-line
            iX[0]=0;
            X=loclx_of_coord(iX);
            su3_put_to_id(Pline[X]); //if glb_t=0 -> P-line=Identity
            if(glb_coord_of_loclx[X][0]>0 && glb_coord_of_loclx[X][0]<=(glb_size[0]/2-1))
	      { //if 0<t<=T/2-1 
		unsafe_su3_prod_su3(aux,Pline[X],plf_dw_buffer);
		su3_copy(Pline[X],aux);
	      }
            for(iX[0]=1;iX[0]<T;iX[0]++)
	      {
                int t=iX[0];
                X=loclx_of_coord(iX);
                //This is a "local" P-line
                if(glb_coord_of_loclx[X][0]<=(glb_size[0]/2-1))
		  {
		    su3_copy(Pline[X],plf[t]);
		    if(rank_coord[0]>0)
		      {
			unsafe_su3_prod_su3(aux,Pline[X],plf_dw_buffer);
			su3_copy(Pline[X],aux);
		      }
		  }
	      }

	    //backward P-line
	   for(iX[0]=T-1;iX[0]>=0;iX[0]--)
	     {
	       int t=iX[0];
	       X=loclx_of_coord(iX);
	       if(glb_coord_of_loclx[X][0]>(glb_size[0]/2-1))
		 {
		   su3_copy(Pline[X],plb[t+1]);
		   if(rank_coord[0]<(nrank_dir[0]-1))
		     {
		       unsafe_su3_prod_su3(aux,Pline[X],plb_up_buffer);
		       su3_copy(Pline[X],aux);
		     }
		 }
	     }	
	}
  
  nissa_free(U0);nissa_free(plf);nissa_free(plb);
}

void Pline_forward(su3 *Pline, quad_su3 *conf)
{
  communicate_lx_quad_su3_borders(conf);
  
  int X,iX[4],T=loc_size[0];
  
  su3 *U0=nissa_malloc("U0",loc_size[0]+1,su3);
  su3 *plf=nissa_malloc("plf",loc_size[0]+1,su3);
  
  MPI_Status status_f;
  int tag_f=0;
  su3 aux,plf_dw_buffer;
  
  //local loop
  for(iX[1]=0;iX[1]<loc_size[1];iX[1]++)
    for(iX[2]=0;iX[2]<loc_size[2];iX[2]++)
      for(iX[3]=0;iX[3]<loc_size[3];iX[3]++)
	{
	  for(iX[0]=0;iX[0]<T;iX[0]++)
	    {
	      X=loclx_of_coord(iX);
	      unsafe_su3_hermitian(U0[iX[0]+1],conf[X][0]);
	    }
           su3_copy(plf[1],U0[1]);
           for(int t=2;t<=T;t++) unsafe_su3_prod_su3(plf[t],U0[t],plf[t-1]);

           if(rank_coord[0]<nrank_dir[0]-1)
	     {
	       su3_copy(plf_dw_buffer,plf[T]);
	       MPI_Send(plf_dw_buffer,1,MPI_SU3,rank_neighup[0],tag_f,MPI_COMM_WORLD);
	     }
           if(rank_coord[0]>0) MPI_Recv(plf_dw_buffer,1,MPI_SU3,rank_neighdw[0],tag_f,MPI_COMM_WORLD,&status_f);

	   iX[0]=0;
	   X=loclx_of_coord(iX);
	   su3_put_to_id(Pline[X]);
	   if(rank_coord[0]>0)
	     {
	       unsafe_su3_prod_su3(aux,Pline[X],plf_dw_buffer);
	       su3_copy(Pline[X],aux);
	     }

           for(iX[0]=1;iX[0]<T;iX[0]++)
	     {
	       X=loclx_of_coord(iX);
	       //This is a "local" P-line
	       su3_copy(Pline[X],plf[iX[0]]);
	       //Here I acummulate the global P-line 
	       if(rank_coord[0]>0)
		 {
		   unsafe_su3_prod_su3(aux,Pline[X],plf_dw_buffer);
		   su3_copy(Pline[X],aux);
		 }
	     }
        }

  nissa_free(U0);nissa_free(plf);
}

void Pline_backward(su3 *Pline, quad_su3 *conf)
{
  communicate_lx_quad_su3_borders(conf);
  
  int X,iX[4],T=loc_size[0];
  
  su3 *U0=nissa_malloc("U0",loc_size[0]+1,su3);
  su3 *plb=nissa_malloc("plb",loc_size[0]+1,su3);
  
  MPI_Status status_b;
  int tag_b=1;
  su3 aux,plb_up_buffer;
  
  //local loop
  for(iX[1]=0;iX[1]<loc_size[1];iX[1]++)
    for(iX[2]=0;iX[2]<loc_size[2];iX[2]++)
      for(iX[3]=0;iX[3]<loc_size[3];iX[3]++)
	{
	  for(iX[0]=0;iX[0]<T;iX[0]++)
	    {
	      X=loclx_of_coord(iX);
	      su3_copy(U0[iX[0]+1],conf[X][0]);
            }
            su3_copy(plb[T],U0[T]);
            for(int t=1;t<T;t++) unsafe_su3_prod_su3(plb[T-t],U0[T-t],plb[T-t+1]);

           if(rank_coord[0]<nrank_dir[0]-1)
	     MPI_Recv(plb_up_buffer,1,MPI_SU3,rank_neighup[0],tag_b,MPI_COMM_WORLD,&status_b);
	   
           if(rank_coord[0]>0)
	     {
	       su3_copy(plb_up_buffer,plb[1]);
	       MPI_Send(plb_up_buffer,1,MPI_SU3,rank_neighdw[0],tag_b,MPI_COMM_WORLD);
	     }

            //backward P-line
	   for(iX[0]=T-1;iX[0]>=0;iX[0]--)
	     {
	       int t=iX[0];
	       X=loclx_of_coord(iX);
	       su3_copy(Pline[X],plb[t+1]);
	       if(rank_coord[0]<nrank_dir[0]-1)
		 {
		   unsafe_su3_prod_su3(aux,Pline[X],plb_up_buffer);
		   su3_copy(Pline[X],aux);
		 }
	     }
	}

  nissa_free(U0);nissa_free(plb);
}

//This will calculate 2*a^2*ig*P_{mu,nu}
/*
  ^                   C--<-- B --<--Y 
  |                   |  2  | |  1  | 
  n                   |     | |     | 
  u                   D-->--\X/-->--A 
  |                   D--<--/X\--<--A 
  -----mu---->        |  3  | |  4  | 
  		      |     | |     | 
		      E-->-- F -->--G 
in order to have the anti-simmetric part, use
the routine "Pmunu_term"
*/
void four_leaves(as2t_su3 *Pmunu,quad_su3 *conf)
{
  communicate_lx_quad_su3_edges(conf);
  
  int A,B,C,D,E,F,G;
  int munu;

  su3 temp1,temp2,leaves_summ;

  nissa_loc_vol_loop(X)
    {
      as2t_su3_put_to_zero(Pmunu[X]);

      munu=0;
      for(int mu=0;mu<4;mu++)
	{
	  A=loclx_neighup[X][mu];
	  D=loclx_neighdw[X][mu];
	  
	  for(int nu=mu+1;nu<4;nu++)
	    {
	      B=loclx_neighup[X][nu];
	      F=loclx_neighdw[X][nu];
	      
	      C=loclx_neighup[D][nu];
	      E=loclx_neighdw[D][nu];
	      
	      G=loclx_neighdw[A][nu];
	      
	      //Put to 0 the summ of the leaves
	      su3_put_to_zero(leaves_summ);
	      
	      //Leaf 1
	      unsafe_su3_prod_su3(temp1,conf[X][mu],conf[A][nu]);         //    B--<--Y 
	      unsafe_su3_prod_su3_dag(temp2,temp1,conf[B][mu]);           //    |  1  | 
	      unsafe_su3_prod_su3_dag(temp1,temp2,conf[X][nu]);           //    |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);	                  //    X-->--A 
	      
	      //Leaf 2
	      unsafe_su3_prod_su3_dag(temp1,conf[X][nu],conf[C][mu]);     //    C--<--B
	      unsafe_su3_prod_su3_dag(temp2,temp1,conf[D][nu]);           //    |  2  | 
	      unsafe_su3_prod_su3(temp1,temp2,conf[D][mu]);		  //    |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);		          //    D-->--X
	      
	      //Leaf 3
	      unsafe_su3_dag_prod_su3_dag(temp1,conf[D][mu],conf[E][nu]);  //   D--<--X
	      unsafe_su3_prod_su3(temp2,temp1,conf[E][mu]);		   //   |  3  | 
	      unsafe_su3_prod_su3(temp1,temp2,conf[F][nu]);		   //   |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);		           //   E-->--F
	      
	      //Leaf 4
	      unsafe_su3_dag_prod_su3(temp1,conf[F][nu],conf[F][mu]);       //  X--<--A 
	      unsafe_su3_prod_su3(temp2,temp1,conf[G][nu]);                 //  |  4  | 
	      unsafe_su3_prod_su3_dag(temp1,temp2,conf[X][mu]);             //  |     |  
	      su3_summ(leaves_summ,leaves_summ,temp1);                      //  F-->--G 
	      
	      su3_copy(Pmunu[X][munu],leaves_summ);
	      
	      munu++;
	    }
	}
    }
  
  set_borders_invalid(Pmunu);
}

//takes the anti-simmetric part of the four-leaves
void Pmunu_term(as2t_su3 *Pmunu,quad_su3 *conf)
{
  four_leaves(Pmunu,conf);
  
  //calculate U-U^dagger
  nissa_loc_vol_loop(X)
    for(int munu=0;munu<6;munu++)
      {
	//bufferized antisimmetrization
	su3 leaves_summ;
	memcpy(leaves_summ,Pmunu[X][munu],sizeof(su3));
	
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++)
	    {
	      Pmunu[X][munu][ic1][ic2][0]=(leaves_summ[ic1][ic2][0]-leaves_summ[ic2][ic1][0])/4;
	      Pmunu[X][munu][ic1][ic2][1]=(leaves_summ[ic1][ic2][1]+leaves_summ[ic2][ic1][1])/4;
	    }
      }
}

//apply the chromo operator to the passed spinor site by site (not yet fully optimized)
void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,as2t_su3 Pmunu,spincolor in)
{
  color temp_d1;
  
  for(int d1=0;d1<4;d1++)
    {
      color_put_to_zero(out[d1]);
      for(int imunu=0;imunu<6;imunu++)
	{
	  unsafe_su3_prod_color(temp_d1,Pmunu[imunu],in[smunu_pos[d1][imunu]]);
	  for(int c=0;c<3;c++) complex_summ_the_prod(out[d1][c],smunu_entr[d1][imunu],temp_d1[c]);
	}
    }
}

//apply the chromo operator to the passed spinor to the whole volume
void unsafe_apply_chromo_operator_to_spincolor(spincolor *out,as2t_su3 *Pmunu,spincolor *in)
{
  nissa_loc_vol_loop(ivol) unsafe_apply_point_chromo_operator_to_spincolor(out[ivol],Pmunu[ivol],in[ivol]);
}

//apply the chromo operator to the passed colorspinspin
//normalization as in ape next
void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin *out,as2t_su3 *Pmunu,colorspinspin *in)
{
  spincolor temp1,temp2;
  
  nissa_loc_vol_loop(ivol)
    {
      //Loop over the four source dirac indexes
      for(int id_source=0;id_source<4;id_source++) //dirac index of source
	{
	  //Switch the color_spinspin into the spincolor.
	  get_spincolor_from_colorspinspin(temp1,in[ivol],id_source);
	  
	  unsafe_apply_point_chromo_operator_to_spincolor(temp2,Pmunu[ivol],temp1);
	  
	  //Switch back the spincolor into the colorspinspin
	  put_spincolor_into_colorspinspin(out[ivol],temp2,id_source);
	}
    }
  
  //invalidate borders
  set_borders_invalid(out);
}

//apply the chromo operator to the passed su3spinspin
//normalization as in ape next
void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin *out,as2t_su3 *Pmunu,su3spinspin *in)
{
  spincolor temp1,temp2;
  
  nissa_loc_vol_loop(ivol)
    {
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
    }

  //invalidate borders
  set_borders_invalid(out);
}

//measure the topological charge site by site
void local_topological_charge(double *charge,quad_su3 *conf)
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
  for(int iperm=0;iperm<3;iperm++)
    {
      //take the index of the two plans
      int ip0=plan_id[iperm][0];
      int ip1=plan_id[iperm][1];
      
      nissa_loc_vol_loop(ivol)
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

//average the topological charge
double average_topological_charge(quad_su3 *conf)
{
  double *charge=nissa_malloc("charge",loc_vol,double);
  
  //compute local charge
  local_topological_charge(charge,conf);
  
  //average over local volume
  double ave_charge=0;
  nissa_loc_vol_loop(ivol)
    ave_charge+=charge[ivol];
  
  nissa_free(charge);
  
  //return the reduction over all nodes
  return glb_reduce_double(ave_charge);
}

//wrapper for eos case
double average_topological_charge(quad_su3 **eo_conf)
{
  //convert to lex
  quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  paste_eo_parts_into_lx_conf(lx_conf,eo_conf);
  
  double charge=average_topological_charge(lx_conf);
  
  nissa_free(lx_conf);
  
  return charge;
}
