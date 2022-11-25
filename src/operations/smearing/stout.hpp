#ifndef _STOUT_HPP
#define _STOUT_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <sstream>

#include "base/field.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //structure to store data for stouting
  struct stout_link_staples
  {
    su3 C;
    su3 Omega;
    su3 Q;
  };
  
  //parameters to stout
  struct stout_pars_t
  {
    int nlevels;
    double rho;
    
    int def_nlevels(){return 0;}
    double def_rho(){return 0;}
    
    int master_fprintf(FILE *fout,int full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(int full=false)
    {
      std::ostringstream os;
      if(full or is_nonstandard())
	{
	  os<<"Stout\n";
	  if(full or nlevels!=def_nlevels()) os<<" NLevels\t=\t"<<nlevels<<"\n";
	  if(full or rho!=def_rho()) os<<" Rho\t\t=\t"<<rho<<"\n";
	}
      return os.str();
    }
    
    int is_nonstandard()
    {
      return
	nlevels!=def_nlevels() or
	rho!=def_rho();
    }
    
    stout_pars_t() :
      nlevels(def_nlevels()),
      rho(def_rho()) {}
  };
  
  //eo
  void stout_smear_whole_stack(eo_ptr<quad_su3> *out,eo_ptr<quad_su3> in,stout_pars_t *stout_pars,const which_dir_t& dirs=all_dirs);
  
  void stout_smear(EoField<quad_su3>& ext_out,
		   const EoField<quad_su3>& ext_in,
		   const stout_pars_t& stout_pars,
		   const which_dir_t& dirs=all_dirs);
  
  void stout_smear_single_level(eo_ptr<quad_su3> out,eo_ptr<quad_su3> ext_in,double rho,const which_dir_t& dirs=all_dirs);
  CUDA_HOST_AND_DEVICE void stout_smear_compute_staples(stout_link_staples *out,eo_ptr<quad_su3> conf,int p,int A,int mu,double rho);
  CUDA_HOST_AND_DEVICE void stout_smear_compute_weighted_staples(su3 staples,eo_ptr<quad_su3> conf,int p,int A,int mu,double rho);
  void stout_smear_conf_stack_allocate(eo_ptr<quad_su3> **out,eo_ptr<quad_su3> in,int nlev);
  void stout_smear_conf_stack_free(eo_ptr<quad_su3> **out,int nlev);
  void stouted_force_remap(eo_ptr<quad_su3> F,eo_ptr<quad_su3> *sme_conf,stout_pars_t *stout_pars);
  void stouted_force_remap_step(eo_ptr<quad_su3> *F,eo_ptr<quad_su3> *conf,double rho);
  //lx
  
  void stout_smear_whole_stack(std::vector<LxField<quad_su3>*>& out,
			       const LxField<quad_su3>& in,
			       const stout_pars_t* stout_pars,
			       const which_dir_t& dirs=all_dirs);
  
  void stout_smear(quad_su3 *ext_out,quad_su3 *ext_in,stout_pars_t *stout_pars,const which_dir_t& dirs=all_dirs);
  void stout_smear_single_level(quad_su3 *out,quad_su3 *ext_in,double rho,const which_dir_t& dirs=all_dirs);
  CUDA_HOST_AND_DEVICE void stout_smear_compute_staples(stout_link_staples *out,quad_su3 *conf,int p,int A,int mu,double rho);

  CUDA_HOST_AND_DEVICE void stout_smear_compute_staples(stout_link_staples *out,
							const LxField<quad_su3>& conf,
							const int& A,
							const int& mu,
							const double& rho);
  
  void stout_smear_conf_stack_allocate(std::vector<LxField<quad_su3>*>& out,
				       LxField<quad_su3>& in,
				       const int& nlev);
  
  void stout_smear_conf_stack_free(std::vector<LxField<quad_su3>*>& out);
  
  void stouted_force_remap(quad_su3 *F,quad_su3 **sme_conf,stout_pars_t *stout_pars);
  
  void stouted_force_remap(LxField<quad_su3>& F,
			   const std::vector<LxField<quad_su3>*>& sme_conf,
			   const stout_pars_t* stout_pars);
  
  /////////////////////////////////////////////////////////////////
  
  /// Compute the lambda entering the force remapping
  template <typename T0,
	    typename T1,
	    typename T2>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void stouted_force_compute_Lambda(const T0& Lambda,
				    const T1& U,
				    const T2& F,
				    hermitian_exp_ingredients *ing)
  {
    complex b[2][3];
    
    if(ing->c1<4e-3)
      {
	double c0=ing->c0,c1=ing->c1;
	b[1][0][RE]=-c0/360;
	b[1][0][IM]=-1.0/6*(1-c1/20*(1-c1/42));
	b[0][0][RE]=0;
	b[0][0][IM]=c0/120*(1-c1/21);
	b[1][1][RE]=1.0/24*(1-c1/15*(1-3*c1/112));
	b[1][1][IM]=-c0/2520;
	b[0][1][RE]=-c0/360*(1-3*c1/56);
	b[0][1][IM]=-1.0/6*(1-c1/10*(1.0-c1/28));
	b[1][2][RE]=0.5*c0/10080;
	b[1][2][IM]=0.5*(1.0/60*(1-c1/21*(1-c1/48)));
	b[0][2][RE]=0.5*(1.0/12*(1-2*c1/30*(1-3*c1/112))); 
	b[0][2][IM]=0.5*(-c0/1260*(1-c1/24));
      }
    else
      {
	//copy back the stored variables
	const double u=ing->u;
	const double w=ing->w;
	const double xi0w=ing->xi0w;
	const double cu=ing->cu,c2u=ing->c2u;
	const double su=ing->su,s2u=ing->s2u;
	const double cw=ing->cw;
	
	//compute additional variables
	const double u2=u*u,w2=w*w;
	double xi1w; //eq. (67)
	if(fabs(w)<0.05) xi1w=-(1-w2*(1-w2*(1-w2/54)/28)/10)/3;
	else xi1w=cw/w2-sin(w)/(w2*w);
	
	//eqs. (60-65)
	complex r[2][3]=
	  {{{2*c2u*u+s2u*(-2*u2+2*w2)+2*cu*u*(8*cw+3*u2*xi0w+w2*xi0w)+su*(-8*cw*u2+18*u2*xi0w+2*w2*xi0w),
	     -8*cw*(2*su*u+cu*u2)+2*(s2u*u+c2u*u2-c2u*w2)+2*(9*cu*u2-3*su*u*u2+cu*w2 -su*u*w2)*xi0w},
	    {2*c2u-4*s2u*u+su*(2*cw*u+6*u*xi0w)+cu*(-2*cw+3*u2*xi0w-w2*xi0w),
	     2*s2u+4*c2u*u+2*cw*(su+cu*u)+(6*cu*u-3*su*u2+su*w2)*xi0w},
	    {-2*s2u+cw*su-3*(su+cu*u)*xi0w,
	     2*c2u+cu*cw+(-3*cu+3*su*u)*xi0w}},
	   {{-2*c2u+2*cw*su*u+2*su*u*xi0w-8*cu*u2*xi0w+6*su*u*u2*xi1w,
	     2*(-s2u+4*su*u2*xi0w+cu*u*(cw+xi0w+3*u2*xi1w))},
	    {2*cu*u*xi0w+su*(-cw-xi0w+3*u2*xi1w),
	     -2*su*u*xi0w-cu*(cw+xi0w-3*u2*xi1w)},
	    {cu*xi0w-3*su*u*xi1w,
	     -(su*xi0w)-3*cu*u*xi1w}}};
	
	//change sign to the f if needed, eq. (34) - changed again below
	if(ing->sign!=0)
	  {
	    ing->f[0][IM]*=-1;
	    ing->f[1][RE]*=-1;
	    ing->f[2][IM]*=-1;
	  }
	
	//compute b
	double t1=9*u2-w2,t2=1/(2*t1*t1);
	for(int j=0;j<3;j++)
	  {
	    //eq. (57-58)
	    for(int ri=0;ri<2;ri++)
	      {
		b[0][j][ri]=(2*u*r[0][j][ri]+(3*u2-w2)*r[1][j][ri]-2*(15*u2+w2)*ing->f[j][ri])*t2;
		b[1][j][ri]=(r[0][j][ri]-3*u*r[1][j][ri]-24*u*ing->f[j][ri])*t2;
	      }
	    
	    //take into account the sign of c0, eq. (70)
	    if(ing->sign!=0)
	      for(int i=0;i<2;i++)
		{
		  //change the sign to real or imag part
		  int ri=(i+j+1)%2;
		  b[i][j][ri]=-b[i][j][ri];
		}
	  }
	
	//change back sign to the f if needed, eq. (34)
	if(ing->sign!=0)
	  {
	    ing->f[0][IM]*=-1;
	    ing->f[1][RE]*=-1;
	    ing->f[2][IM]*=-1;
	  }
      }
    
    //compute B eq. (69)
    su3 B[2];
    for(int i=0;i<2;i++)
      {
	su3_put_to_diag(B[i],b[i][0]);
	su3_summ_the_prod_complex(B[i],ing->Q, b[i][1]);
	su3_summ_the_prod_complex(B[i],ing->Q2,b[i][2]);
      }
    
    //compute Gamma, eq. (74)
    su3 Gamma;
    //compute U*Sigma', to be used many times
    su3 aux;
    unsafe_su3_prod_su3(aux,U,F);
    //compute the trace of U*Sigma'*B[j]
    complex we[2];
    for(int j=0;j<2;j++) trace_su3_prod_su3(we[j],aux,B[j]);
    //first term
    unsafe_su3_prod_complex(Gamma,ing->Q,we[0]);
    //first term
    su3_summ_the_prod_complex(Gamma,ing->Q2,we[1]);
    //third term
    su3_summ_the_prod_complex(Gamma,aux,ing->f[1]);
    //fourth and fifth term
    su3 temp;
    unsafe_su3_prod_su3  (temp,ing->Q,aux);
    su3_summ_the_prod_su3(temp,aux,ing->Q);
    su3_summ_the_prod_complex(Gamma,temp,ing->f[2]);
    
    //compute Lambda eq. (73)
    unsafe_su3_traceless_hermitian_part(Lambda,Gamma);
  }
  
}

#endif
