#include "nissa.h"

#define rad2 1.414213562373095048801688724209

int *send_of_iter[2],*recv_of_iter[2];

void reorder_gaugeconf(su3 *reord,quad_su3 *orig)
{
  int j=0,k=2*loc_vol;
  for(int i=0;i<loc_volr;i++)
    {
      su3_copy(reord[j++],orig[i][0]);
      su3_copy(reord[j++],orig[i][1]);
  
      su3_copy(reord[k++],orig[i][2]);
      su3_copy(reord[k++],orig[i][3]);
    }
}

void setup_senders_receivers_and_gauge(int *send_of_iter[2],su3 *reord_U,quad_su3 *U)
{
  for(int eo_or_oe=0;eo_or_oe<2;eo_or_oe++) send_of_iter[eo_or_oe]=malloc(4*loc_volr*sizeof(int));
  
  int nbord_dir[4]={0,0,0,0};
  for(int loceo=0;loceo<loceo;loceo++)
    {
      for(int idir=0;idir<4;idir++)
	{
	  //Application of forward projector and derivative will be required from backward link, in even position
	  send_of_iter[0][loceo]=8*loce_neighdw[loceo][idir]+2*idir;
	  send_of_iter[1][loceo]=8*loco_neighdw[loceo][idir]+2*idir;

	  memcpy(reord_U[0][0][4*loceo+idir],U[loclx_of_loco[loce_neighdw[loceo][idir]]][idir],sizeof(su3));
	  memcpy(reord_U[1][0][4*loceo+idir],U[loclx_of_loce[loceo]][idir],sizeof(su3));
	  
	  memcpy(reord_U[0][1][4*loceo+idir],U[loclx_of_loce[loco_neighdw[loceo][idir]]][idir],sizeof(su3));
	  memcpy(reord_U[1][1][4*loceo+idir],U[loclx_of_loco[loceo]][idir],sizeof(su3));
	  
	  //Application of backward projector will be required from forward link, in odd position
	  send_of_iter[0][loceo]=8*loce_neighup[loceo][idir]+2*idir+1;
	  send_of_iter[1][loceo]=8*loco_neighup[loceo][idir]+2*idir+1;
	}
    }
}

void apply_DOE_or_DEO(spincolor *out,spincolor *in,su3 *reord_U[2][2],int eo_or_oe)
{
  color psi,psi2;
  spincolor rs;
  spincolor *s=in;
  redspincolor *phi=allocate_redspincolor(loc_vol/2+bord_vol,"phi");
  su3 *U;
  
  int *loc_send_of_iter=send_of_iter[eo_or_oe];
  
  U=reord_U[0][eo_or_oe];
  
  int ix=0;
  for(int i=0;i<loc_volr;i++)
    {
      spincolor_copy(rs,*s);
      s++;
      
      /*********************** direction +0 ************************/
      
      int loc_send=(*(loc_send_of_iter++));

      color_summ(psi,rs[0],rs[2]);
      color_summ(psi2,rs[1],rs[3]);

      unsafe_su3_prod_color(phi[loc_send][0],*U,psi);
      unsafe_su3_prod_color(phi[loc_send][1],*U,psi2);
            
      U++;ix++;
    
      /*********************** direction -0 ************************/
      
      loc_send=(*(loc_send_of_iter++));
      
      color_subt(phi[loc_send][0],rs[0],rs[2]);
      color_subt(phi[loc_send][1],rs[1],rs[3]);

      ix++;

      /*********************** direction +1 ************************/
      
      loc_send=(*(loc_send_of_iter++));
      
      color_isumm(psi,rs[0],rs[3]);
      color_isumm(psi2,rs[1],rs[2]);
      
      unsafe_su3_prod_color(phi[loc_send][0],*U,psi);
      unsafe_su3_prod_color(phi[loc_send][1],*U,psi2);
      
      U++;ix++;

      /*********************** direction -1 ************************/
      
      loc_send=(*(loc_send_of_iter++));
      
      color_isubt(phi[loc_send][0],rs[0],rs[3]);
      color_isubt(phi[loc_send][1],rs[1],rs[2]);

      ix++;

      /*********************** direction +2 ************************/

      loc_send=(*(loc_send_of_iter++));
      
      color_summ(psi,rs[0],rs[3]);
      color_subt(psi2,rs[1],rs[2]);

      unsafe_su3_prod_color(phi[loc_send][0],*U,psi);
      unsafe_su3_prod_color(phi[loc_send][1],*U,psi2);
      
      U++;ix++;

      /*********************** direction -2 ************************/
      
      loc_send=(*(loc_send_of_iter++));
 
      color_subt(phi[loc_send][0],rs[0],rs[3]);
      color_summ(phi[loc_send][1],rs[1],rs[2]);

      ix++;

      /*********************** direction +3 ************************/
      
      loc_send=(*(loc_send_of_iter++));
      
      color_isumm(psi,rs[0],rs[2]);
      color_isubt(psi2,rs[1],rs[3]);      

      unsafe_su3_prod_color(phi[loc_send][0],*U,psi);
      unsafe_su3_prod_color(phi[loc_send][1],*U,psi2);

      U++;ix++;

      /*********************** direction -3 ************************/
      
      loc_send=(*(loc_send_of_iter++));
      
      color_isubt(phi[loc_send][0],rs[0],rs[2]);
      color_isumm(phi[loc_send][1],rs[1],rs[3]);
      
      ix++;
      
      /************************ end of loop ************************/
    }
  
  s=out;

  U=reord_U[1][eo_or_oe];

  ix=0;

  int loc_recv=0;
  for(int i=0;i<loc_volr;i++)
    {
      /*********************** direction +0 ************************/
	
      color_copy(rs[0],phi[loc_recv][0]);
      color_copy(rs[2],phi[loc_recv][0]);
      color_copy(rs[1],phi[loc_recv][1]);
      color_copy(rs[3],phi[loc_recv][1]);
      
      loc_recv++;ix++;
	
      /*********************** direction -0 ************************/
	
      unsafe_su3_dag_prod_color(psi,*U,phi[loc_recv][0]);
      unsafe_su3_dag_prod_color(psi2,*U,phi[loc_recv][1]);
	
      summassign_color(rs[0],psi);
      subtassign_color(rs[2],psi);
      summassign_color(rs[1],psi2);
      subtassign_color(rs[3],psi2);
	
      loc_recv++;ix++;U++;
	
      /*********************** direction +1 ************************/
      
      summassign_color(rs[0],phi[loc_recv][0]);
      subtassign_icolor(rs[3],phi[loc_recv][0]);
      summassign_color(rs[1],phi[loc_recv][1]);
      subtassign_icolor(rs[2],phi[loc_recv][1]);
      
      loc_recv++;ix++;
	
      /*********************** direction -1 ************************/
      
      unsafe_su3_dag_prod_color(psi, *U,phi[loc_recv][0]);
      unsafe_su3_dag_prod_color(psi2,*U,phi[loc_recv][1]);
      
      summassign_color(rs[0],psi);
      summassign_icolor(rs[3],psi);
      summassign_color(rs[1],psi2);
      summassign_icolor(rs[2],psi2);
      
      loc_recv++;U++;ix++;
	
      /*********************** direction +2 ************************/
      
      summassign_color(rs[0],phi[loc_recv][0]);
      summassign_color(rs[3],phi[loc_recv][0]);
      
      summassign_color(rs[1],phi[loc_recv][1]);
      subtassign_color(rs[2],phi[loc_recv][1]);
      
      loc_recv++;ix++;
      
      /*********************** direction -2 ************************/
      
      unsafe_su3_dag_prod_color(psi,*U,phi[loc_recv][0]);
      unsafe_su3_dag_prod_color(psi2,*U,phi[loc_recv][1]);

      summassign_color(rs[0],psi);
      subtassign_color(rs[3],psi);
      summassign_color(rs[1],psi2);
      summassign_color(rs[2],psi2);
      
      loc_recv++;U++;ix++;
	
      /*********************** direction +3 ************************/
      
      summassign_color(rs[0],phi[loc_recv][0]);
      subtassign_icolor(rs[2],phi[loc_recv][0]);
      
      summassign_color(rs[1],phi[loc_recv][1]);
      summassign_icolor(rs[3],phi[loc_recv][1]);
      
      loc_recv++;ix++;
      
      /*********************** direction -3 ************************/
      
      unsafe_su3_dag_prod_color(psi,*U,phi[loc_recv][0]);
      unsafe_su3_dag_prod_color(psi2,*U,phi[loc_recv][1]);
      
      color_summ((*s)[0],rs[0],psi);
      color_isumm((*s)[2],rs[2],psi);
      color_summ((*s)[1],rs[1],psi2);
      color_isubt((*s)[3],rs[3],psi2);
      
      loc_recv++;U++;ix++;s++;
    }
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();

  if(narg<3 && rank==0)
    {
      fprintf(stderr,"Use: %s L T\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  int L=atoi(arg[1]);
  int T=atoi(arg[2]);

  //Init the MPI grid 
  init_grid(L,T);
  init_random(0);

  ///////////////////////////////////////////

  set_eo_geometry();
  
  int *send_of_iter[2];
  quad_su3 *conf=allocate_quad_su3(loc_vol+bord_vol,"conf");
  su3 *reord_U[2][2];
  for(int eo_or_oe=0;eo_or_oe<2;eo_or_oe++)
    {
      send_of_iter[eo_or_oe]=(int*)malloc(sizeof(int)*loc_volr);
      for(int OT=0;OT<2;OT++) reord_U[OT][eo_or_oe]=allocate_su3(loc_vol*8,"reord_U");
    }
  
  setup_senders_receivers_and_guauge(send_of_iter,reord_U,conf);
  
  for(int i=0;i<loc_vol*2;i++)
    if(rank==0)
      printf("%d %d\n",i,send_of_iter[0][i]);
  
  //////////////////////////////////////////

  close_nissa();

  return 0;
}
