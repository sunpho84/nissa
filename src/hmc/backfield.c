#pragma once

//initialize an u(1) field to unity
void init_backfield_to_id(quad_u1 **S)
{
  for(int par=0;par<2;par++)
    nissa_loc_volh_loop(ivol)
      for(int mu=0;mu<4;mu++)
	{
	  S[par][ivol][mu][0]=1;
	  S[par][ivol][mu][1]=0;
	}
}

/*
//add the staggered phases to a background field
void add_stagphases_to_backfield(quad_u1 **S)
{
  for(int ilx=0;ilx<loc_vol;ilx++)
    {
      int d[4];
      d[1]=0;
      d[2]=d[1]+glb_coord_of_loclx[ilx][1];
      d[3]=d[2]+glb_coord_of_loclx[ilx][2];
      d[0]=d[3]+glb_coord_of_loclx[ilx][3];
      
      int par=loclx_parity[ilx];
      int ieo=loceo_of_loclx[ilx];
      
      //define each dir
      for(int mu=0;mu<4;mu++)
	for(int ri=0;ri<2;ri++)
	  S[par][ieo][mu][ri]*=(d[mu]%2==1) ? -1 : 1;
    }
}

//add temporary anti-periodic boundary condition to a background field
void add_antiperiodic_bc_to_backfield(quad_u1 **S)
{
  for(int par=0;par<2;par++)
    for(int ieo=0;ieo<loc_volh;ieo++)
      if(glb_coord_of_loclx[loclx_of_loceo[par][ieo]][0]==glb_size[0]-1)
	for(int ri=0;ri<2;ri++)
	  S[par][ieo][0][ri]*=-1;
}
*/

//multpiply the configuration for an additional u(1) field
void add_backfield_to_conf(quad_su3 **conf,quad_u1 **u1)
{
  master_printf("Adding backfield\n");
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  safe_su3_prod_complex(conf[par][ivol][mu],conf[par][ivol][mu],u1[par][ivol][mu]);
      set_borders_invalid(conf[par]);
    }
}

//multpiply the configuration for an the conjugate of an u(1) field
void rem_backfield_from_conf(quad_su3 **conf,quad_u1 **u1)
{
  master_printf("Removing backfield\n");
  for(int par=0;par<2;par++)
    {
      nissa_loc_volh_loop(ivol)
	for(int mu=0;mu<4;mu++)
	  safe_su3_prod_conj_complex(conf[par][ivol][mu],conf[par][ivol][mu],u1[par][ivol][mu]);
      set_borders_invalid(conf[par]);
    }
}
