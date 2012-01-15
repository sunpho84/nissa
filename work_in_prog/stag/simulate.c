#include "nissa.h"
#include "new_struct.c"

double beta=5.38;
double am=0.025;

int L,T;
quad_su3 *eo_conf[2];
color *pf[2];
int npf=2;
quad_su3 *H[2],*F[2];
quad_u1 *u1b[2];

//allocate a new rat approx
rat_approx* rat_approx_malloc(int nterms,char *name)
{
  rat_approx *out=malloc(sizeof(rat_approx));
  memcpy(out->name,name,20);
  out->nterms=nterms;
  out->poles=malloc(sizeof(double)*2*nterms);
  out->weights=out->poles+nterms;
  
  return out;
}

//free a rational approx
void rat_approx_free(rat_approx *ptr)
{
  free(ptr->poles);
  free(ptr);
}

//print a rational approximation
void master_print_rat_approx(rat_approx *ptr)
{
  master_printf("Rational approximation %s:\n",ptr->name);
  master_printf("  const: %lg\n",ptr->con);
  master_printf("  nterms: %d\n",ptr->nterms);
  for(int iterm=0;iterm<nterms;iterm++)
    master_printf("   %d) pole: %lg, coef: %lg\n",iterm,ptr->poles[iterm],ptr->coef[iterm]);
}

//initialize an u(1) field to unity
void init_backfield_to_id(quad_u1 **S)
{
  for(int par=0;par<2;par++)
    for(int ieo=0;ieo<loc_vol/2;ieo++)
      for(int mu=0;mu<4;mu++)
	{
	  S[par][ieo][mu][0]=1;
	  S[par][ieo][mu][1]=0;
	}
}

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
    for(int ieo=0;ieo<loc_vol/2;ieo++)
      if(glb_coord_of_loclx[loclx_of_loceo[par][ieo]][0]==glb_size[0]-1)
	for(int ri=0;ri<2;ri++)
	  S[par][ieo][0][ri]*=-1;
}

//multpiply the configuration for an additional u(1) field
void add_backfield_to_conf(quad_su3 **conf,quad_u1 **u1)
{
  for(int par=0;par<2;par++)
    for(int ieo=0;ieo<loc_vol/2;ieo++)
      for(int mu=0;mu<4;mu++)
	safe_su3_prod_complex(conf[par][ieo][mu],conf[par][ieo][mu],u1[par][ieo][mu]);
}

//multpiply the configuration for an the conjugate of an u(1) field
void rem_backfield_to_conf(quad_su3 **conf,quad_u1 **u1)
{
  for(int par=0;par<2;par++)
    for(int ieo=0;ieo<2;ieo++)
      for(int mu=0;mu<4;mu++)
	safe_su3_prod_conj_complex(conf[par][ieo][mu],conf[par][ieo][mu],u1[par][ieo][mu]);
}

//initialize the simulation
void init_simulation()
{
  //basic mpi initialization
  init_nissa();

  //set sizes
  L=16;
  T=4;
  
  //Init the MPI grid 
  init_grid(T,L);
  
  //initialize the local random generators
  start_loc_rnd_gen(0);
  
  //allocate the conf
  eo_conf[0]=nissa_malloc("conf",loc_vol+loc_bord,quad_su3);
  eo_conf[1]=eo_conf[0]+(loc_vol+loc_bord)/2;
  
  //allocate the u1 background field
  u1b[0]=nissa_malloc("u1_back",loc_vol+loc_bord,quad_u1);
  u1b[1]=u1b[0]+(loc_vol+loc_bord)/2;
  
  //allocate pseudo-fermions
  for(int iq=0;iq<2;iq++)
    pf[iq]=nissa_malloc("pf",(loc_vol+loc_bord)/2,color);
  
  //allocate the momenta
  H[0]=nissa_malloc("H",loc_vol+loc_bord,quad_su3);
  H[1]=H[0]+(loc_vol+loc_bord)/2;
  
  //allocate the force
  F[0]=nissa_malloc("F",loc_vol+loc_bord,quad_su3);
  F[1]=F[0]+(loc_vol+loc_bord)/2;
  
  //read the conf
  read_ildg_conf_and_split_into_eo_parts(eo_conf,"conf");
  
  //initialize background field to stagg phases and anti-periodic bc
  init_backfield_to_id(u1b);
  add_stagphases_to_backfield(u1b);
  add_antiperiodic_bc_to_backfield(u1b);
  
  //allocate rational approximation for pseudo-fermions generation
  ra_pf=rat_approx_malloc(15,"pforig");
}

//generate momenta using guassian hermitean matrix generator
void generate_momenta()
{
  for(int eo=0;eo<2;eo++)
    for(int ivol=0;ivol<loc_vol/2;ivol++)
      for(int mu=0;mu<4;mu++)
	herm_put_to_gauss(H[eo][ivol][mu],&(loc_rnd_gen[ivol]),1);
}

//debug
void read_e_color(color *e,char *path)
{
  color *lx=nissa_malloc("lx",loc_vol,color);
  
  read_color(lx,path);
  take_e_part_of_lx_color(e,lx,loc_vol);
  
  nissa_free(lx);
}
void read_o_color(color *o,char *path)
{
  color *lx=nissa_malloc("lx",loc_vol,color);
  
  read_color(lx,path);
  take_o_part_of_lx_color(o,lx,loc_vol);
  
  nissa_free(lx);
}

//generate pseudo-fermions using color vector generator
void generate_pseudo_fermions()
{
  /*
  color *temp=nissa_malloc("temp",loc_vol/2,color);  
  //generate the random field
  for(int ivol=0;ivol<loc_vol/2;ivol++)
    color_put_to_gauss(temp[ivol],&(loc_rnd_gen[ivol]),1);
  nissa_free(temp);
  */
  
  color *read_rnd=nissa_malloc("rnd",(loc_vol+loc_bord)/2,color);
  color *read_phi=nissa_malloc("read_phi",(loc_vol+loc_bord)/2,color);
  color *comp_phi=nissa_malloc("comp_phi",(loc_vol+loc_bord)/2,color);
  
  add_backfield_to_conf(eo_conf,u1b);
  communicate_eo_gauge_borders(eo_conf[0],eo_conf[1]);
  
  read_e_color(read_rnd,"rnd");
  communicate_ev_color_borders(read_rnd);
  
  
  double con=1.9135313739066195;
  double coef[15]={-1.99205872380355886E-007,-8.17161635223895208E-007,-2.64017630352392784E-006,-8.20070280923918561E-006,-2.52434088690784991E-005,-7.75084159688460014E-005,-2.37848908569994309E-004,-7.30230885251800148E-004,-2.24648131931426366E-003,-6.95100598099882352E-003,-2.18488934822891161E-002,-7.16925229561548721E-002,-0.26496783456064565,-1.3840888450676481,-24.990790193256750};
  
  double mass2[15]={3.16155247975239949E-006,1.89211624893373918E-005,6.44801500601987701E-005,1.88879830821950712E-004,5.26032179337476789E-004,1.43893286604875850E-003,3.91088319579291008E-003,1.06077212226758488E-002,2.87754464106065588E-002,7.82470690950233805E-002,0.21433502691613024,0.59921220654378748,1.7746138289708084,6.2621765269663721,44.035100776270461};
  
  for(int imass=0;imass<15;imass++) mass2[imass]+=0.025*0.025;
  summ_src_and_all_inv_stD2ee_cgmm2s(comp_phi,read_rnd,eo_conf,con,mass2,coef,15,1000000,1.e-10,1.e-10,0);
  
  //apply_stD2ee(comp_phi,eo_conf,read_phi,0.025,read_rnd);
  
  read_e_color(read_phi,"phie");
  
  double diff=0;
  for(int ivol=0;ivol<loc_vol/2;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	{
	  master_printf("%lg %lg\n",comp_phi[ivol][ic][ri],read_phi[ivol][ic][ri]);
	  diff+=pow(comp_phi[ivol][ic][ri]-read_phi[ivol][ic][ri],2);
	}
  diff/=loc_vol/2*2*3;

  master_printf("%lg\n",sqrt(diff));
  
  nissa_free(comp_phi);
  nissa_free(read_phi);
  nissa_free(read_rnd);
}

//finalize everything
void close_simulation()
{
  nissa_free(eo_conf[0]);
  nissa_free(H[0]);
  nissa_free(F[0]);
  nissa_free(u1b[0]);
  for(int iq=0;iq<2;iq++)
    nissa_free(pf[iq]);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_simulation();
  
  ///////////////////////////////////////
  
  generate_momenta();
  generate_pseudo_fermions();
  
  ///////////////////////////////////////
  
  close_simulation();
  
  return 0;
}
