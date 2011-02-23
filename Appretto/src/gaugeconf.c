#pragma once

//rotate the gauge configuration anti-clockwise by 90 degrees
//on one of the spatial dirs

/*

   .---.---.---.        .---.---.---.     
   |   2       |        |           |            	
   .   A1  .   .        .   .   .   .     
   |           |        |       1   |     	
   .   .   .   .        .   C2' B   .     
   |           |        |           |     	
   O---.---.---.        O---.---.---.     

   d2
   O d1

*/

void ac_rotate_gauge_conf(quad_su3 *out,quad_su3 *in,int axis)
{
  int d1=1+(axis-1+1)%3;
  int d2=1+(axis-1+2)%3;
  
  if(rank==0) printf("Direction to rotate: %d %d\n",d1,d2);

  if(rank==0 && rank_tot>1)
    {
      fprintf(stderr,"Rotation works only on single node!\n");
      MPI_Abort(MPI_COMM_WORLD,0);
    }
  
  if(rank==0 && glb_size[d1]!=glb_size[d2])
    {
      fprintf(stderr,"Rotation works only if dir %d and %d have the same size!\n",glb_size[d1],glb_size[d2]);
      MPI_Abort(MPI_COMM_WORLD,0
);
    }
  
  int xA[4],xB[4],xC[4];
  for(xA[0]=0;xA[0]<glb_size[0];xA[0]++)
    for(xA[1]=0;xA[1]<glb_size[1];xA[1]++)
      for(xA[2]=0;xA[2]<glb_size[2];xA[2]++)
	for(xA[3]=0;xA[3]<glb_size[3];xA[3]++)
	  {
	    xC[0]=xB[0]=xA[0];
	    xC[axis]=xB[axis]=xA[axis];

	    xC[d2]=xB[d2]=xA[d1];
	    xC[d1]=xB[d1]=glb_size[d2]-1-xA[d2];
	    
	    xC[d1]=(xC[d1]-1+glb_size[d1])%glb_size[d1];
	    
	    int A=glblx_of_coord(xA);
	    int B=glblx_of_coord(xB);
	    int C=glblx_of_coord(xC);

	    su3_copy(out[B][0],in[A][0]);
	    su3_copy(out[B][axis],in[A][axis]);
	    su3_copy(out[B][d2],in[A][d1]);
	    for(int ic1=0;ic1<3;ic1++)
	      for(int ic2=0;ic2<3;ic2++)
		{
		  out[C][d1][ic1][ic2][0]=+in[A][d2][ic2][ic1][0];
		  out[C][d1][ic1][ic2][1]=-in[A][d2][ic2][ic1][1];
		}
	  }
}

//shift the gauge configuration
void shift_gauge_conf_down(quad_su3 *conf,int *amount)
{
  quad_su3 *temp_conf=allocate_quad_su3(loc_vol+loc_bord,"temp_conf");
  quad_su3 *supp_conf=allocate_quad_su3(loc_vol,"supp_conf");

  //initial copy inside temp_conf
  memcpy(temp_conf,conf,sizeof(quad_su3)*loc_vol);

  for(int idir=0;idir<4;idir++)
    for(int ix=0;ix<amount[idir];ix++)
      {
        communicate_gauge_borders(temp_conf);
        
        for(int ivol=0;ivol<loc_vol;ivol++)
          {
            int orig=loclx_neighup[ivol][idir];
            quad_su3_copy(supp_conf[ivol],temp_conf[orig]);
          }

        memcpy(temp_conf,supp_conf,sizeof(quad_su3)*loc_vol);
      }

  memcpy(conf,temp_conf,sizeof(quad_su3)*loc_vol);

  free(temp_conf);free(supp_conf);
}

//put boundary conditions on the gauge conf
void put_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges)
{
  complex theta[4];
  for(int idir=0;idir<4;idir++)
    {
      theta[idir][0]=cos(theta_in_pi[idir]*PI/glb_size[idir]);
      theta[idir][1]=sin(theta_in_pi[idir]*PI/glb_size[idir]);
    }
  
  int nsite=loc_vol;
  if(putonbords) nsite+=loc_bord;
  if(putonedges) nsite+=loc_edge;

  for(int loc_site=0;loc_site<nsite;loc_site++)
    for(int idir=0;idir<4;idir++) safe_su3_prod_complex(conf[loc_site][idir],conf[loc_site][idir],theta[idir]);
}

void rem_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges)
{
  double meno_theta_in_pi[4]={-theta_in_pi[0],-theta_in_pi[1],-theta_in_pi[2],-theta_in_pi[3]};
  put_boundaries_conditions(conf,meno_theta_in_pi,putonbords,putonedges);
}

//Adapt the border condition
void adapt_theta(quad_su3 *conf,double *old_theta,double *put_theta,int putonbords,int putonedges)
{
  double diff_theta[4];
  int adapt=0;

  for(int idir=0;idir<4;idir++)
    {
      adapt=adapt || (old_theta[idir]!=put_theta[idir]);
      diff_theta[idir]=put_theta[idir]-old_theta[idir];
      old_theta[idir]=put_theta[idir];
    }
  
  if(adapt)
    {
      if(rank==0) printf("Necesarry to add boundary condition: %f %f %f %f\n",diff_theta[0],diff_theta[1],diff_theta[2],diff_theta[3]);
      put_boundaries_conditions(conf,diff_theta,putonbords,putonedges);
    }
}
