#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <fcntl.h>
#include <unistd.h>

#include <base/field.hpp>
#include <geometry/geometry_eo.hpp>
#include <geometry/geometry_lx.hpp>
#include <new_types/su3_op.hpp>
#include <routines/mpi_routines.hpp>
#include <threads/threads.hpp>

#ifndef EXTERN_RANDOM
# define EXTERN_RANDOM extern
#endif

//random number generator table length
#define RAN2_NTAB 32

namespace nissa
{
  //Random types
  const int nrnd_type=7;
  enum rnd_t{RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_UNIF,RND_Z2,RND_Z3,RND_Z4,RND_GAUSS};
  const char rnd_t_str[nrnd_type][20]={"AllPlusOne","AllMinusOne","Unif","Z2","Z3","Z4","Gauss"};
  
  //Source type
  enum source_t{POINT_SOURCE,UNDILUTED_SOURCE,COLOR_DILUTED_SOURCE,SPIN_DILUTED_SOURCE,SPINCOLOR_DILUTED_SOURCE};
  
  //The structure for the random generator
  struct rnd_gen
  {
    int idum;
    int idum2;
    int iv[RAN2_NTAB];
    int iy;
  };
  
  //random generator stuff
  EXTERN_RANDOM rnd_gen glb_rnd_gen;
  EXTERN_RANDOM bool glb_rnd_gen_inited;
  CUDA_MANAGED EXTERN_RANDOM rnd_gen *loc_rnd_gen;
  EXTERN_RANDOM bool loc_rnd_gen_inited;
  
  rnd_t convert_str_to_rnd_t(const char *str);
  void color_put_to_gauss(color H,rnd_gen *gen,double sigma);
  void convert_text_to_rnd_gen(rnd_gen *gen,const char *text);
  void convert_rnd_gen_to_text(char *text,rnd_gen *gen,int size);
  CUDA_HOST_AND_DEVICE double rnd_get_unif(rnd_gen *gen,double min,double max);
  CUDA_HOST_AND_DEVICE int rnd_get_pm_one(rnd_gen *gen);
  
  //generate a spindiluted vector according to the passed type
  void generate_colorspindiluted_source(LxField<su3spinspin>& source,
					const rnd_t& rtype,
					const int& twall);
  
  //generate a spindiluted vector according to the passed type
  void generate_spindiluted_source(LxField<colorspinspin>& source,
				   const rnd_t& rtype,
				   const int& twall);
  
  //generate an undiluted vector according to the passed type
  void generate_undiluted_source(LxField<spincolor>& source,
				 const rnd_t& rtype,
				 const int& twall);
  
  //generate a fully undiluted source
  void generate_fully_undiluted_lx_source(LxField<color>& source,
					  const rnd_t& rtype,
					  const int& twall,
					  const int& dir=0);
  
  //eo version
  void generate_fully_undiluted_eo_source(EvenOrOddField<color>& source,
					  const rnd_t& rtype,
					  const int& twall,
					  const int& par,
					  const int& dir=0);
  
  void generate_fully_undiluted_eo_source(EoField<color>& source,
					  const rnd_t& rtype,
					  const int& twall,
					  const int& dir=0);
  
  //same for spincolor
  void generate_fully_undiluted_eo_source(EvenOrOddField<spincolor>& source,
					  const rnd_t& rtype,
					  const int& twall,
					  const int& par,
					  const int& dir=0);
  
  void generate_fully_undiluted_eo_source(EoField<spincolor>& source,
					  const rnd_t& rtype,
					  const int& twall,
					  const int& dir=0);
  
  //generate a delta source
  void generate_delta_source(LxField<su3spinspin>& source,
			     const coords_t& x);
  
  void generate_delta_eo_source(EoField<su3>& source,
				const coords_t& x);

  // void generate_delta_eo_source(eo_ptr<su3> source,int *x);
  // void generate_delta_source(su3spinspin *source,int *x);
  // void generate_colorspindiluted_source(su3spinspin *source,enum rnd_t rtype,int twall);
  // inline void generate_spincolordiluted_source(su3spinspin *source,enum rnd_t rtype,int twall)
  // {generate_colorspindiluted_source(source,rtype,twall);}
  // void generate_spindiluted_source(colorspinspin *source,enum rnd_t rtype,int twall);
  // void generate_undiluted_source(spincolor *source,enum rnd_t rtype,int twall);
  // void generate_fully_undiluted_lx_source(color *source,enum rnd_t rtype,int twall,int dir=0);
  
  // void generate_fully_undiluted_eo_source(EvenOrOddField<color>& source,
  // 					  const rnd_t& rtype,
  // 					  const int& twall,
  // 					  const int& par,
  // 					  const int& dir=0);
  
  // void generate_fully_undiluted_eo_source(EoField<color>& source,
  // 					  const rnd_t& rtype,
  // 					  const int& twall,
  // 					  const int& par,
  // 					  const int& dir=0);
  
  // void generate_fully_undiluted_eo_source(EvenOrOddField<spincolor>& source,
  // 					  const rnd_t& rtype,
  // 					  const int& twall,
  // 					  const int& par,
  // 					  const int& dir=0);
  
  // void generate_fully_undiluted_eo_source(eo_ptr<spincolor> source,enum rnd_t rtype,int twall,int dir=0);
  
  CUDA_HOST_AND_DEVICE void herm_put_to_gauss(su3& H,rnd_gen *gen,double sigma);
  void rnd_fill_pm_one_loc_vector(double *v,int nps);
  void rnd_fill_unif_loc_vector(double *v,int dps,double min,double max);
  coords_t generate_random_coord();
  double rnd_get_gauss_double(rnd_gen *gen,double ave=0,double sig=1);
  
  //return a Z2 complex
  template <typename C>
  CUDA_HOST_AND_DEVICE void rnd_get_Z2(C&& out,
				       rnd_gen *gen)
  {
    out[0]=rnd_get_pm_one(gen);
    out[1]=0;
  }
  
  template <typename C>
  CUDA_HOST_AND_DEVICE inline void rnd_get_Z3(C&& out,
					      rnd_gen *gen)
  {
    rnd_get_ZN(out,gen,3);
  }
  
  //return a Z4 complex
  template <typename C>
  CUDA_HOST_AND_DEVICE void rnd_get_Z4(C&& out,
				       rnd_gen *gen)
  {
    out[0]=rnd_get_pm_one(gen)/(double)RAD2;
    out[1]=rnd_get_pm_one(gen)/(double)RAD2;
  }
  
  //return a ZN complex
  template <typename C>
  CUDA_HOST_AND_DEVICE void rnd_get_ZN(C&& out,
				       rnd_gen *gen,
				       const int& N)
  {
    complex_iexp(out,2*M_PI*(int)rnd_get_unif(gen,0,N)/N);
  }
  
  //return a gaussian complex with sigma=sig/sqrt(2)
  template <typename C>
  CUDA_HOST_AND_DEVICE void rnd_get_gauss_complex(C&& out,
						  rnd_gen *gen,
						  const complex& ave,
						  const double& sig)
  {
    const double one_by_sqrt2=0.707106781186547;
    double norm=sig*one_by_sqrt2;
    double q,r;
    
    r=sqrt(-2*log(1-rnd_get_unif(gen,0,1)));
    q=2*M_PI*rnd_get_unif(gen,0,1);
    
    out[0]=r*cos(q)*norm+ave[0];
    out[1]=r*sin(q)*norm+ave[1];
  }
  
  template <typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void comp_get_rnd(C&& out,
		    rnd_gen *gen,
		    const enum rnd_t& rtype)
  {
    switch(rtype)
      {
      case RND_ALL_PLUS_ONE: complex_put_to_real(out,+1);                   break;
      case RND_ALL_MINUS_ONE:complex_put_to_real(out,-1);                   break;
      case RND_UNIF:         complex_put_to_real(out,rnd_get_unif(gen,0,1));break;
      case RND_Z2:           rnd_get_Z2(out,gen);                           break;
      case RND_Z3:           rnd_get_Z3(out,gen);                           break;
      case RND_Z4:           rnd_get_Z4(out,gen);                           break;
      case RND_GAUSS:        rnd_get_gauss_complex(out,gen,{0.0,0.0},1);    break;
      }
  }
  
  void start_glb_rnd_gen(const char *text);
  void start_glb_rnd_gen(int seed);
  void start_loc_rnd_gen(int seed);
  void start_loc_rnd_gen(const char *mess);
  void start_rnd_gen(rnd_gen *out,int seed);
  void stop_loc_rnd_gen();
  void su3_find_heatbath(su3 out,su3 in,su3 staple,double beta,int nhb_hits,rnd_gen *gen);
  
    /// Put a matrix to random used passed random generator
  template <typename U>
  CUDA_HOST_AND_DEVICE
  void su3_put_to_rnd(U&& u_ran,
		      rnd_gen &rnd)
  {
    su3_put_to_id(u_ran);
    
    for(size_t i1=0;i1<NCOL;i1++)
      for(size_t i2=i1+1;i2<NCOL;i2++)
	{
	  //generate u0,u1,u2,u3 random on the four dim. sphere
	  const double u0=rnd_get_unif(&rnd,-1,1);
	  const double alpha=sqrt(1-u0*u0);
	  const double phi=rnd_get_unif(&rnd,0,2*M_PI);
	  const double costheta=rnd_get_unif(&rnd,-1,1);
	  const double sintheta=sqrt(1-costheta*costheta);
	  const double u3=alpha*costheta;
	  const double u1=alpha*sintheta*cos(phi);
	  const double u2=alpha*sintheta*sin(phi);
	  
	  //define u_l as unit matrix ...
	  su3 u_l;
	  su3_put_to_id(u_l);
	  
	  //... and then modify the elements in the chosen su(2) subgroup
	  u_l[i1][i1][RE]=u0;
	  u_l[i1][i1][IM]=u3;
	  u_l[i1][i2][RE]=u2;
	  u_l[i1][i2][IM]=u1;
	  u_l[i2][i1][RE]=-u2;
	  u_l[i2][i1][IM]=u1;
	  u_l[i2][i2][RE]=u0;
	  u_l[i2][i2][IM]=-u3;
	  
	  safe_su3_prod_su3(u_ran,u_l,u_ran);
	}
  }
}

#undef EXTERN_RANDOM

#endif
