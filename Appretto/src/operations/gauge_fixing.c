#pragma once

//determine the gauge transformation bringing to temporal gauge with T-1 timeslice diferent from id
void find_temporal_gauge_fixing_matr(su3 *fixm,quad_su3 *u)
{
  int loc_slice_area=loc_size[1]*loc_size[2]*loc_size[3];
  su3 *buf=NULL;
  
  //if the number of processors in the 0 dir is greater than 1 allocate room for border
  if(nproc_dir[0]>1) buf=appretto_malloc("buf",loc_slice_area,su3);

  //if we are on first proc slice put to identity the t=0 slice, otherwise receive it from previous proc slice
  if(proc_coord[0]==0)
    {
      for(int ivol=0;ivol<loc_vol;ivol++)
	if(glb_coord_of_loclx[ivol][0]==0)
	  su3_put_to_id(fixm[ivol]);
    }
  else
    if(nproc_dir[0]>1)
      MPI_Recv((void*)fixm,loc_slice_area,MPI_SU3,rank_neighdw[0],252,cart_comm,MPI_STATUS_IGNORE);
  
  //now go ahead along t
  int c[4];
  //loop over spatial slice
  for(c[1]=0;c[1]<loc_size[1];c[1]++)
    for(c[2]=0;c[2]<loc_size[2];c[2]++)
      for(c[3]=0;c[3]<loc_size[3];c[3]++)
	{
	  //bulk
	  for(c[0]=1;c[0]<loc_size[0];c[0]++)
	    {
	      int icurr=loclx_of_coord(c);
	      c[0]--;int iback=loclx_of_coord(c);c[0]++;
	      
	      su3_prod_su3(fixm[icurr],fixm[iback],u[iback][0]);
	    }
	  //border
	  if(nproc_dir[0]>1)
	    {
	      c[0]=loc_size[0]-1;int iback=loclx_of_coord(c);
	      c[0]=0;int icurr=loclx_of_coord(c);
	      
	      su3_prod_su3(buf[icurr],fixm[iback],u[iback][0]);
	    }
	  
	}
  
  //if we are not on last slice of processor send g to next slice
  if(proc_coord[0]!=(nproc_dir[0]-1) && nproc_dir[0]>1)
    MPI_Send((void*)buf,loc_slice_area,MPI_SU3,rank_neighup[0],252,cart_comm);

  if(nproc_dir[0]>1) appretto_free(buf);
}

//apply a gauge transformation to the conf
void gauge_transform_conf(quad_su3 *uout,su3 *g,quad_su3 *uin)
{
  if(rank_tot>1)
    {
      communicate_su3_borders(g);
      communicate_gauge_borders(uin);
    }
  
  su3 temp;
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int mu=0;mu<4;mu++)
      {
	su3_prod_su3_dag(temp,uin[ivol][mu],g[loclx_neighup[ivol][mu]]);
	su3_prod_su3(uout[ivol][mu],g[ivol],temp);
      }
}
