#pragma once

//compute the action of the momenta
double momenta_action(quad_su3 **H)
{
  double loc_action=0;
  
  //summ the square of H
  for(int par=0;par<2;par++)
    nissa_loc_volh_loop(ivol)
      for(int mu=0;mu<4;mu++)
	loc_action+=real_part_of_trace_su3_prod_su3_dag(H[par][ivol][mu],H[par][ivol][mu]);
  
  //global reducton
  double glb_action;
  MPI_Allreduce(&loc_action,&glb_action,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  return glb_action*0.5;  
}
