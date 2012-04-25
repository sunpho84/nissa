#include <math.h>

#include "../new_types/complex.h"
#include "../new_types/new_types_definitions.h"
#include "../base/debug.h"
#include "../base/global_variables.h"
#include "../base/routines.h"
#include "../base/vectors.h"
#include "../geometry/geometry_lx.h"

#include "random.h"

//Initialize a random number generator
void start_rnd_gen(rnd_gen *out,int seed)
{
  const int im1=2147483563,ia1=40014;
  const int iq1=53668,ir1=12211;
  int j,k;
  
  //initialization
  out->idum=seed;
  out->idum=max_int(out->idum+1,1);
  out->idum2=out->idum;
  for(j=ran2_ntab+7;j>=0;j--)
    {
      k=out->idum/iq1;
      out->idum=ia1*(out->idum-k*iq1)-k*ir1;
      if(out->idum<0) out->idum+=im1;
          if(j<ran2_ntab) out->iv[j]=out->idum;
    }
  out->iy=out->iv[0];
}

//initialize the global random generator
void start_glb_rnd_gen(int seed)
{
  if(nissa_glb_rnd_gen_inited==1) crash("global random generator already initialized");
  start_rnd_gen(&(glb_rnd_gen),seed);
  
  nissa_glb_rnd_gen_inited=1;
  master_printf("Global random generators initialized with seed: %d\n",seed);
}

//Initialize the grid of local random number generator
void start_loc_rnd_gen(int seed)
{
  if(nissa_loc_rnd_gen_inited==1) crash("local random generator already initialized!");
  
  //check the grid to be initiaized
  if(loc_vol==0) crash("grid not initalized!");
  
  //Generate the true seed
  if(nissa_glb_rnd_gen_inited==0) start_glb_rnd_gen(seed);
  int internal_seed=(int)rnd_get_unif(&glb_rnd_gen,0,RAND_MAX);
  
  //allocate the grid of random generator, one for site
  loc_rnd_gen=nissa_malloc("Loc_rnd_gen",loc_vol,rnd_gen);
  for(int ivol=0;ivol<loc_vol;ivol++) start_rnd_gen(&(loc_rnd_gen[ivol]),internal_seed+glblx_of_loclx[ivol]);
  
  nissa_loc_rnd_gen_inited=1;
  master_printf("Grid of local random generators initialized with internal seed: %d\n",internal_seed);
}

//Stop grid of local random generators
void stop_loc_rnd_gen()
{
  if(nissa_loc_rnd_gen_inited==0) crash("local random generator not initialized");
  master_printf("Stopping local random generators\n");
  
  nissa_free(loc_rnd_gen);
  nissa_loc_rnd_gen_inited=0;
}

//standard ran2 from numerical recipes
double rnd_get_unif(rnd_gen *gen,double min,double max)
{
  const int im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014,ia2=40692;
  const int iq1=53668,iq2=52774,ir1=12211,ir2=3791,ndiv=1+imm1/ran2_ntab;
  const double am=1.0/im1,eps=1.2e-7,rnmx=1-eps;
  int j,k;
  double out;

  k=gen->idum/iq1;
  gen->idum=ia1*(gen->idum-k*iq1)-k*ir1;
  if(gen->idum<0) gen->idum+=im1;
  
  k=gen->idum2/iq2;
  gen->idum2=ia2*(gen->idum2-k*iq2)-k*ir2;
  if(gen->idum2<0) gen->idum2+=im2;
  
  j=gen->iy/ndiv;
  gen->iy=gen->iv[j]-gen->idum2;
  gen->iv[j]=gen->idum;
  if(gen->iy<0) gen->iy+=imm1;
  
  out=min_double(am*gen->iy,rnmx);
  
  return out*(max-min)+min;
}

//return a numer between 0 and 1
int rnd_get_pm_one(rnd_gen *gen)
{
  double r=rnd_get_unif(gen,0,1);
  if(r>0.5) return 1;
  else return -1;
}

//return a Z2 complex
void rnd_get_Z2(complex out,rnd_gen *gen)
{
  out[0]=rnd_get_pm_one(gen);
  out[1]=0;
}

//return a Z4 complex
void rnd_get_Z4(complex out,rnd_gen *gen)
{
  out[0]=rnd_get_pm_one(gen)/RAD2;
  out[1]=rnd_get_pm_one(gen)/RAD2;
}

//return a gaussian real
double rnd_get_gauss(rnd_gen *gen,double ave,double sig)
{
  double q,r,x;
  static double y;
  static int flag=1;
  
  if(flag)
    {
      r=sqrt(-2*log(1-rnd_get_unif(gen,0,1)));
      q=2*M_PI*rnd_get_unif(gen,0,1);
      
      x=r*cos(q);
      y=r*sin(q);
    }
  else x=y;
  
  flag=!flag;
  
  return x*sig+ave;
}

//return a gaussian complex whose module has sigma sig
void rnd_get_gauss_complex(complex out,rnd_gen *gen,complex ave,double sig)
{
  const double one_by_sqrt2=0.707106781186547;
  out[0]=rnd_get_gauss(gen,ave[0],sig*one_by_sqrt2);
  out[1]=rnd_get_gauss(gen,ave[1],sig*one_by_sqrt2);
}

//return a complex number of appropriate type
void comp_get_rnd(complex out,rnd_gen *gen,enum rnd_type rtype)
{
  complex z={0,0};
  switch(rtype)
    {
    case RND_ALL_PLUS_ONE:  out[0]=1;                     out[1]=0;break;
    case RND_ALL_MINUS_ONE: out[0]=-1;                    out[1]=0;break;
    case RND_UNIF:          out[0]=rnd_get_unif(gen,0,1); out[1]=0;break;
    case RND_Z2:            rnd_get_Z2(out,gen);                   break;
    case RND_Z4:            rnd_get_Z4(out,gen);                   break;
    case RND_GAUSS:         rnd_get_gauss_complex(out,gen,z,1);    break;
    }
}

//fill a grid of vectors with numbers between 0 and 1
void rnd_fill_unif_loc_vector(double *v,int dps,double min,double max)
{
  nissa_loc_vol_loop(ivol)
    for(int i=0;i<dps;i++)
      v[ivol*dps+i]=rnd_get_unif(&(loc_rnd_gen[ivol]),min,max);
}

//return a grid of +-x numbers
void rnd_fill_pm_one_loc_vector(double *v,int nps)
{
  nissa_loc_vol_loop(ivol)
    for(int i=0;i<nps;i++)
      v[ivol*nps+i]=rnd_get_pm_one(&(loc_rnd_gen[ivol]));
  
  set_borders_invalid(v);
}

//generate a spindiluted vector according to the passed type
void generate_spindiluted_source(colorspinspin *source,enum rnd_type rtype,int twall)
{ //reset
  memset(source,0,sizeof(colorspinspin)*loc_vol);
  
  //compute normalization norm: spat vol if twall>=0, glb_vol else
  int norm2=glb_vol;
  if(twall>=0) norm2/=glb_size[0];
  double inv_sqrt_norm2=1.0/sqrt(norm2);
  
  nissa_loc_vol_loop(ivol)
    if(glb_coord_of_loclx[ivol][0]==twall||twall<0)
      for(int ic=0;ic<3;ic++)
	{ //generate the id_sink==id_source==0 entry
	  comp_get_rnd(source[ivol][ic][0][0],&(loc_rnd_gen[ivol]),rtype);
	  complex_prod_double(source[ivol][ic][0][0],source[ivol][ic][0][0],inv_sqrt_norm2);
	  for(int d=1;d<4;d++) //copy the other three dirac indexes
	    memcpy(source[ivol][ic][d][d],source[ivol][ic][0][0],sizeof(complex));
	}
  
  set_borders_invalid(source);
}

//generate an undiluted vector according to the passed type
void generate_undiluted_source(spincolor *source,enum rnd_type rtype,int twall)
{ //reset
  memset(source,0,sizeof(spincolor)*loc_vol);
  
  nissa_loc_vol_loop(ivol)
    if(glb_coord_of_loclx[ivol][0]==twall||twall<0)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  comp_get_rnd(source[ivol][id][ic],&(loc_rnd_gen[ivol]),rtype);
  
  set_borders_invalid(source);
}

//generate a delta source
void generate_delta_source(su3spinspin *source,int *x)
{ //reset
  memset(source,0,sizeof(su3spinspin)*loc_vol);

  int islocal=1,lx[4];
  for(int idir=0;idir<4;idir++)
    {
      lx[idir]=x[idir]-rank_coord[idir]*loc_size[idir];
      islocal&=(lx[idir]>=0);
      islocal&=(lx[idir]<loc_size[idir]);
    }
  
  if(islocal)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
        source[loclx_of_coord(lx)][ic][ic][id][id][0]=1;
  
  set_borders_invalid(source);
}
