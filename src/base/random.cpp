#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>

#define EXTERN_RANDOM
# include "base/random.hpp"

#include "base/debug.hpp"
#include "base/field.hpp"
#include "base/random.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"

namespace nissa
{
  // template <typename IMin,
  // 	    typename IMax,
  // 	    typename F>
  // __global__
  // void cudageneric_kernel(const IMin &min,
  // 			   const IMax &max,
  // 			   F f)
  // {
  //   const auto i=min+blockIdx.x*blockDim.x+threadIdx.x;
  //   if(i<max)
  //     f(loclx_neighdw[0][0]);
  // }
  
  // void rrr()
  // {
  //   const dim3 block_dimension(NUM_THREADS);
  //   const dim3 grid_dimension((1+block_dimension.x-1)/block_dimension.x);
  //   cudageneric_kernel<<<grid_dimension,block_dimension>>>(0, 1, [=] __host__ __device__(int){});
  // }
  
  CUDA_HOST_AND_DEVICE double rnd_get_unif(rnd_gen *gen,double min,double max);
  
  //initialize a random number generator
  void start_rnd_gen(rnd_gen *out,int seed)
  {
    const int im1=2147483563,ia1=40014;
    const int iq1=53668,ir1=12211;
    int j,k;
    
    //initialization
    out->idum=seed;
    out->idum=std::max(out->idum+1,1);
    out->idum2=out->idum;
    for(j=RAN2_NTAB+7;j>=0;j--)
      {
	k=out->idum/iq1;
	out->idum=ia1*(out->idum-k*iq1)-k*ir1;
	if(out->idum<0) out->idum+=im1;
	if(j<RAN2_NTAB) out->iv[j]=out->idum;
      }
    out->iy=out->iv[0];
  }
  
  //print all the entries of the random generator into a string
  void convert_rnd_gen_to_text(char *text,rnd_gen *gen,int size)
  {
    int *ptr=(int*)gen;
    memset(text,0,size);
    
    //init output text
    snprintf(text,size,"%d",ptr[0]);
    
    //append all the elements
    for(int i=1;i<RAN2_NTAB+3;i++)
      {
	char temp[20];
	snprintf(temp,20," %d",ptr[i]);
	strncat(text,temp,size);
      }
    if(strlen(text)==1024) crash("use larger text");
  }
  
  //read all the entries of the random generator from a string
  void convert_text_to_rnd_gen(rnd_gen *gen,const char *text)
  {
    int *ptr=(int*)gen;
    
    for(int i=0;i<RAN2_NTAB+3;i++)
      {
	char temp[20];
	if(sscanf(text,"%s",temp)!=1) crash("while reading element %d from %s",i,text);
	text+=1+strlen(temp);
	if(sscanf(temp,"%d",ptr+i)!=1) crash("while converting to int %s",temp);
      }
  }
  
  //initialize the global random generator
  void start_glb_rnd_gen(int seed)
  {
    if(glb_rnd_gen_inited==1) crash("global random generator already initialized");
    start_rnd_gen(&(glb_rnd_gen),seed);
    
    glb_rnd_gen_inited=1;
    master_printf("Global random generators initialized with seed: %d\n",seed);
  }
  
  //init from text
  void start_glb_rnd_gen(const char *text)
  {
    if(glb_rnd_gen_inited==1) crash("global random generator already initialized");
    convert_text_to_rnd_gen(&(glb_rnd_gen),text);
    
    glb_rnd_gen_inited=1;
    master_printf("Global random generators initialized from text\n");
  }
  
  //initialize the grid of local random number generator
  void start_loc_rnd_gen(int seed)
  {
    if(loc_rnd_gen_inited==1) crash("local random generator already initialized!");
    
    //check the grid to be initiaized
    if(locVol==0) crash("grid not initalized!");
    
    //Generate the true seed
    if(glb_rnd_gen_inited==0) start_glb_rnd_gen(seed);
    int internal_seed=(int)rnd_get_unif(&glb_rnd_gen,0,RAND_MAX-glbVol);
    
    //allocate the grid of random generator, one for site
    loc_rnd_gen=nissa_malloc("Loc_rnd_gen",locVol,rnd_gen);
    for(int ivol=0;ivol<locVol;ivol++)
      {
	int loc_seed=internal_seed+glblxOfLoclx[ivol];
	if(loc_seed<0) crash("something went wrong with local seed: %d + %d = %d",internal_seed,glblxOfLoclx[ivol],loc_seed);
	start_rnd_gen(&(loc_rnd_gen[ivol]),loc_seed);
      }
    loc_rnd_gen_inited=1;
    master_printf("Grid of local random generators initialized with internal seed: %d\n",internal_seed);
  }
  
  //init from text
  void start_loc_rnd_gen(const char *text)
  {
    start_glb_rnd_gen(text);
    start_loc_rnd_gen(0); //glb grid already started so does not matter
  }
  
  //stop grid of local random generators
  void stop_loc_rnd_gen()
  {
    if(loc_rnd_gen_inited==0) crash("local random generator not initialized");
    master_printf("Stopping local random generators\n");
    
    nissa_free(loc_rnd_gen);
    loc_rnd_gen_inited=0;
  }
  
  //standard ran2 from numerical recipes
  CUDA_HOST_AND_DEVICE double rnd_get_unif(rnd_gen *gen,double min,double max)
  {
    const int im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014,ia2=40692;
    const int iq1=53668,iq2=52774,ir1=12211,ir2=3791,ndiv=1+imm1/RAN2_NTAB;
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
    
    out=nissa_min(am*gen->iy,rnmx);
    
    return out*(max-min)+min;
  }
  
  //generate a random postion
  coords_t generate_random_coord()
  {
    coords_t c;
    for(int mu=0;mu<NDIM;mu++)
      c[mu]=(int)(rnd_get_unif(&glb_rnd_gen,0,glbSize[mu]));
    
    return c;
  }
  
  //return a numer between 0 and 1
  CUDA_HOST_AND_DEVICE int rnd_get_pm_one(rnd_gen *gen)
  {
    double r=rnd_get_unif(gen,0,1);
    if(r>0.5) return 1;
    else return -1;
  }
  
  //return a gaussian double
  double rnd_get_gauss_double(rnd_gen *gen,double ave,double sig)
  {
    double q,r;
    
    r=sqrt(-2*log(1-rnd_get_unif(gen,0,1)));
    q=2*M_PI*rnd_get_unif(gen,0,1);
    
    return r*cos(q)*sig+ave;
  }
  
  //get the type of random from a string
  rnd_t convert_str_to_rnd_t(const char *str)
  {
    for(int i=0;i<nrnd_type;i++)
      if(strcasecmp(str,rnd_t_str[i])==0)
	return (rnd_t)i;
    
    master_fprintf(stderr,"Error, unknown random string %s, known ones:\n",str);
    for(int i=0;i<nrnd_type;i++)
      master_fprintf(stderr," %s\n",rnd_t_str[i]);
    crash("Choose one of them");
    
    return (rnd_t)RND_ALL_MINUS_ONE;
  }
  
  // //fill a grid of vectors with numbers between 0 and 1
  // void rnd_fill_unif_loc_vector(double* v,int dps,double min,double max)
  // {
  //   NISSA_PARALLEL_LOOP(ivol,0,locVol)
  //     for(int i=0;i<dps;i++)
  // 	v[ivol*dps+i]=rnd_get_unif(&(loc_rnd_gen[ivol]),min,max);
  //   NISSA_PARALLEL_LOOP_END;
    
  //   set_borders_invalid(v);
  // }
  
  //return a grid of +-x numbers
  // void rnd_fill_pm_one_loc_vector(double* v,int nps)
  // {
  //   NISSA_PARALLEL_LOOP(ivol,0,locVol)
  //     for(int i=0;i<nps;i++)
  // 	v[ivol*nps+i]=rnd_get_pm_one(&(loc_rnd_gen[ivol]));
  //   NISSA_PARALLEL_LOOP_END;
    
  //   set_borders_invalid(v);
  // }
  
  //generate a spindiluted vector according to the passed type
  void generate_colorspindiluted_source(LxField<su3spinspin>& source,
					const rnd_t& rtype,
					const int& twall)
  {
    //reset
    source.reset();
    
    PAR(0,locVol,CAPTURE(twall,rtype,TO_WRITE(source)),ivol,
	{
	  if(glbCoordOfLoclx[ivol][0]==twall or twall<0)
	    {
	      comp_get_rnd(source[ivol][0][0][0][0],&(loc_rnd_gen[ivol]),rtype);
	      for(int c=0;c<NCOL;c++)
		for(int d=0;d<NDIRAC;d++)
		  if(c or d)
		    complex_copy(source[ivol][c][c][d][d],source[ivol][0][0][0][0]);
	    }
	});
  }
  
  //generate a spindiluted vector according to the passed type
  void generate_spindiluted_source(LxField<colorspinspin>& source,
				   const rnd_t& rtype,
				   const int& twall)
  {
    source.reset();
    
    PAR(0,locVol,CAPTURE(twall,rtype,TO_WRITE(source)),ivol,
	{
	  if(glbCoordOfLoclx[ivol][0]==twall or twall<0)
	    for(int ic=0;ic<NCOL;ic++)
	      {
		comp_get_rnd(source[ivol][ic][0][0],&(loc_rnd_gen[ivol]),rtype);
		for(int d=1;d<NDIRAC;d++)
		  complex_copy(source[ivol][ic][d][d],source[ivol][ic][0][0]);
	      }
	});
  }
  
  //generate an undiluted vector according to the passed type
  void generate_undiluted_source(LxField<spincolor>& source,
				 const rnd_t& rtype,
				 const int& twall)
  {
    source.reset();
    
    PAR(0,locVol,CAPTURE(twall,rtype,TO_WRITE(source)),ivol,
	{
	  if(glbCoordOfLoclx[ivol][0]==twall or twall<0)
	    for(int id=0;id<NDIRAC;id++)
	      for(int ic=0;ic<NCOL;ic++)
		comp_get_rnd(source[ivol][id][ic],&(loc_rnd_gen[ivol]),rtype);
	});
  }
  
  //generate a fully undiluted source
  void generate_fully_undiluted_lx_source(LxField<color>& source,
					  const rnd_t& rtype,
					  const int& twall,
					  const int& dir)
  {
    source.reset();
    
    PAR(0,locVol,CAPTURE(dir,twall,rtype,TO_WRITE(source)),ilx,
	{
	  if(twall<0 or glbCoordOfLoclx[ilx][dir]==twall)
	    for(int ic=0;ic<NCOL;ic++)
	      comp_get_rnd(source[ilx][ic],&(loc_rnd_gen[ilx]),rtype);
	});
  }
  
  //eo version
  void generate_fully_undiluted_eo_source(EvenOrOddField<color>& source,
					  const rnd_t& rtype,
					  const int& twall,
					  const int& par,
					  const int& dir)
  {
    source.reset();
    
    PAR(0,locVolh,CAPTURE(par,dir,twall,rtype,TO_WRITE(source)),ieo,
	{
	  int ilx=loclx_of_loceo[par][ieo];
	  if(twall<0 or glbCoordOfLoclx[ilx][dir]==twall)
	    for(int ic=0;ic<NCOL;ic++)
	      comp_get_rnd(source[ieo][ic],&(loc_rnd_gen[ilx]),rtype);
	});
  }
  
  void generate_fully_undiluted_eo_source(EoField<color>& source,
					  const rnd_t& rtype,
					  const int& twall,
					  const int& dir)
  {
    for(int par=0;par<2;par++)
      generate_fully_undiluted_eo_source(source[par],rtype,twall,par,dir);
  }
  
  //same for spincolor
  void generate_fully_undiluted_eo_source(EvenOrOddField<spincolor>& source,
					  const rnd_t& rtype,
					  const int& twall,
					  const int& par,
					  const int& dir)
  {
    source.reset();
    
    PAR(0,locVolh,CAPTURE(par,dir,twall,rtype,TO_WRITE(source)),ieo,
	{
	  int ilx=loclx_of_loceo[par][ieo];
	  if(twall<0 or glbCoordOfLoclx[ilx][dir]==twall)
	    for(int id=0;id<NDIRAC;id++)
	      for(int ic=0;ic<NCOL;ic++)
		comp_get_rnd(source[ieo][id][ic],&(loc_rnd_gen[ilx]),rtype);
	});
  }
  
  void generate_fully_undiluted_eo_source(EoField<spincolor>& source,
					  const rnd_t& rtype,
					  const int& twall,
					  const int& dir)
  {
    for(int par=0;par<2;par++)
      generate_fully_undiluted_eo_source(source[par],rtype,twall,par,dir);
  }
  
  //generate a delta source
  void generate_delta_source(LxField<su3spinspin>& source,
			     const coords_t& x)
  {
    //reset
    source.reset();
    
    int islocal=1;
    coords_t lx;
    for(int mu=0;mu<NDIM;mu++)
      {
	lx[mu]=x[mu]-rank_coord[mu]*locSize[mu];
	islocal&=(lx[mu]>=0);
	islocal&=(lx[mu]<locSize[mu]);
      }
    
    if(islocal)
      {
	auto s=source.getWritable();
	
	for(int id=0;id<4;id++)
	  for(int ic=0;ic<NCOL;ic++)
	  s[loclx_of_coord(lx)][ic][ic][id][id][0]=1;
      }
  }
  
  //generate a delta source
  void generate_delta_eo_source(EoField<su3>& source,
				const coords_t& x)
  {
    //reset
    source.reset();
    
    int islocal=1;
    coords_t lx;
    for(int mu=0;mu<NDIM;mu++)
      {
        lx[mu]=x[mu]-rank_coord[mu]*locSize[mu];
        islocal&=(lx[mu]>=0);
        islocal&=(lx[mu]<locSize[mu]);
      }
    
    if(islocal)
      {
	auto s=source.getWritable();
	
	const int ivol=loclx_of_coord(lx);
	su3_put_to_id(s[loclx_parity[ivol]][loceo_of_loclx[ivol]]);
      }
  }
  
  /////////////////////////////// Generate an hermitian matrix ///////////////////////
  
  // A gauss vector has complex components z which are gaussian distributed
  // with <z~ z> = sigma
  void color_put_to_gauss(color& H,rnd_gen *gen,double sigma)
  {
    complex ave={0,0};
    for(size_t ic=0;ic<NCOL;ic++) rnd_get_gauss_complex(H[ic],gen,ave,sigma);
  }
  
  //return a single link after the heatbath procedure
  //routines to be shrunk!
  void su3_find_heatbath(su3 out,su3 in,su3 staple,double beta,int nhb_hits,rnd_gen *gen)
  {
    //compute the original contribution to the action due to the given link 
    su3 prod;
    unsafe_su3_prod_su3_dag(prod,in,staple);
    
    //copy in link to out
    if(out!=in) su3_copy(out,in);
    
    //iterate over heatbath hits
    for(int ihit=0;ihit<nhb_hits;ihit++)
      //scan all the three possible subgroups
      for(int isub_gr=0;isub_gr<3;isub_gr++)
	{
	  //take the part of the su3 matrix
	  double r0,r1,r2,r3;
	  double smod=su2_part_of_su3(r0,r1,r2,r3,prod,isub_gr);
	  
	  //omega is the coefficient of the plaquette, divided by the module of the su2 submatrix for normalization
	  double omega_f=beta/(3*smod);
	  double z_norm=exp(-2*omega_f);
	  omega_f=1/omega_f;
	  
	  double temp_f,z_f,a0;
	  do
	    {
	      double z_temp=(z_norm-1)*rnd_get_unif(gen,0,1)+1;
	      a0     = 1+omega_f*log(z_temp);
	      z_f    = 1-a0*a0;
	      temp_f = sqr(rnd_get_unif(gen,0,1))-z_f;
	    }
	  while(temp_f>0);
	  
	  double x_rat=sqrt(z_f);
	  
	  //generate an su2 matrix
	  double fi=rnd_get_unif(gen,0,2*M_PI);
	  double cteta=rnd_get_unif(gen,-1,1);
	  double steta=sqrt(1-cteta*cteta);
	  
	  double a1=steta*cos(fi)*x_rat;
	  double a2=steta*sin(fi)*x_rat;
	  double a3=cteta*x_rat;
	  
	  double x0 = a0*r0 + a1*r1 + a2*r2 + a3*r3;
	  double x1 = r0*a1 - a0*r1 + a2*r3 - r2*a3;
	  double x2 = r0*a2 - a0*r2 + a3*r1 - r3*a1;
	  double x3 = r0*a3 - a0*r3 + a1*r2 - r1*a2;
	  
	  su2_prodassign_su3(x0,x1,x2,x3,isub_gr,prod);
	  su2_prodassign_su3(x0,x1,x2,x3,isub_gr,out);
	}
  }
}
