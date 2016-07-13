#ifndef _GEOMETRY_VIR_HPP
#define _GEOMETRY_VIR_HPP

#ifndef EXTERN_GEOMETRY_VIR
 #define EXTERN_GEOMETRY_VIR extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

#include <functional>

#include "base/grid.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"

#include "new_types/su3.hpp"
#include "new_types/float_128.hpp"

namespace nissa
{
#if SIMD_TYPE == SIMD_EMU
  const int simd_width=256;
#endif
  const int nvranks_max=simd_width/(8*sizeof(float));
  
  //! switch whether to use vranks or not
  EXTERN_GEOMETRY_VIR int use_vranks;
  //! minimal virtual volume
  EXTERN_GEOMETRY_VIR int vloc_min_vol;
  //! minimal size for the virtual ranks
  EXTERN_GEOMETRY_VIR coords vloc_min_size;
  //! fix the maximal number of virtual ranks
  EXTERN_GEOMETRY_VIR coords fix_nvranks_max;
  
  //! remap to lx from
  void lx_remap_to_virsome_internal(void *out,void *in,int size_per_site,int nel_per_site,int nvranks,int *vrank_of_loclx,int *idx);
  void virsome_remap_to_lx_internal(void *out,void *in,int size_per_site,int nel_per_site,int nvranks,int *vrank_of_loclx,int *vrank_loclx_offset,int *idx);
  
  //! structure to hold a virtual grid
  template<class base_type> struct vranks_grid_t
  {
    //! number of elements inside a simd
    static const int nvranks=simd_width/(8*sizeof(base_type));
    
    //types
    typedef base_type vbase_type[nvranks];
    typedef vbase_type vcomplex[2];
    typedef vcomplex vcolor[NCOL];
    typedef vcolor vsu3[NCOL];
    typedef vcolor vhalfspincolor[NDIRAC/2];
    typedef vhalfspincolor vcolor_halfspincolor[NCOL];
    typedef vcolor_halfspincolor vhalfspincolor_halfspincolor[NDIRAC/2];
    typedef vsu3 vquad_su3[NDIM];
    typedef vsu3 voct_su3[2*NDIM];
    typedef vcolor vspincolor[NDIRAC];
    typedef vcomplex vhalfspin[2];
    typedef vsu3 vclover_term_t[4];
    typedef vhalfspincolor_halfspincolor vinv_clover_term_t[2];
    
    //! initialization flag
    int inited;
    //! hold which directions are parallelized
    coords vparal_dir;
    //! coordinateds of the elements
    coords vrank_coord[nvranks];
    //! number of parallelized directions
    int nvparal_dir;
    //! virtual rank local volume
    int vloc_vol;
    //! volume of the even or odd sublattices
    int vloc_volh;
    //! size of the virtual ranks grid
    coords vloc_size;
    //! number of virtual ranks per direction
    coords nvranks_per_dir;
    //! offset between corresponding elements in vrank 0
    int vrank_loclx_offset[nvranks];
    //! virtual rank hosting a local lattice element
    int *vrank_of_loclx;
    
    //! return the rank of the loceo
    int vrank_of_loceo(int par,int eo)
    {return vrank_of_loclx[loclx_of_loceo[par][eo]];}
    
    //! initialize the grid
    void init()
    {
      if(!inited)
      {
  	inited=true;
	
	//CRASH if eo-geom not inited
	if(!eo_geom_inited)
	  {
	    if(!use_eo_geom) CRASH("eo geometry must be enabled in order to use vir one");
	    CRASH("initialize eo_geometry before vir one");
	  }
	
	//set volume of the virtual node
	vloc_vol=loc_vol/nvranks;
	vloc_volh=vloc_vol/2;
	master_printf("Initializing vgeom for type of %d bits, nvranks=%d\n",(int)sizeof(base_type)*8,nvranks);
	vrank_of_loclx=nissa_malloc("vrank_of_loclx",loc_vol,int);
	
      	//fix the number of minimal node
	partitioning_t vranks_partitioning(loc_vol,nvranks);
	master_printf("Grouping loc_vol %d in %d vranks obtained %u possible combos\n",loc_vol,nvranks,vranks_partitioning.ncombo);
	coords VRPD;
	coords min_loc_size,fix_nvranks;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    min_loc_size[mu]=2;
	    fix_nvranks[mu]=0;
	  }
	master_printf("need to find the optimal here ANNNNNNA\n");
	bool something_found=false;
	while(vranks_partitioning.find_next_valid_partition(VRPD,loc_size,min_loc_size,fix_nvranks))
	  {
	    //set the vir local size
	    //coords VLS;
	    //for(int mu=0;mu<NDIM;mu++) VLS[mu]=nvloc_max_per_dir[mu]/VPD[mu];
	    
	    //if it is the minimal surface (or first valid combo) copy it and compute the border size
	    //if(rel_surf<min_rel_surf||min_rel_surf<0)
	      {
		//min_rel_surf=rel_surf;
		
		for(int mu=0;mu<NDIM;mu++)
		  {
		    nvranks_per_dir[mu]=VRPD[mu];
		    vloc_size[mu]=loc_size[mu]/VRPD[mu];
		  }
		
	master_printf(" proposed local volume\t%d",vloc_size[0]);
	for(int mu=1;mu<NDIM;mu++) master_printf("x%d",vloc_size[mu]);
	master_printf(" = %d\n",vloc_vol);
		something_found=true;
	      }
	  }
	if(!something_found) CRASH("no valid grid partitioning found");
	
	//compute the local virtual and size check whether each direction dir is virtually parallelized
	nvparal_dir=0;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    vloc_size[mu]=loc_size[mu]/nvranks_per_dir[mu];
	    vparal_dir[mu]=(nvranks_per_dir[mu]>1);
	    nvparal_dir+=vparal_dir[mu];
	  }
	
	//create the virtual grid
	for(int i=0;i<nvranks;i++) coord_of_lx(vrank_coord[i],i,nvranks_per_dir);
	
	//print information on virtual local volume
	master_printf("Virtual local volume\t%d",vloc_size[0]);
	for(int mu=1;mu<NDIM;mu++) master_printf("x%d",vloc_size[mu]);
	master_printf(" = %d\n",vloc_vol);
	master_printf("List of virtually parallelized dirs:\t");
	for(int mu=0;mu<NDIM;mu++) if(vparal_dir[mu]) master_printf("%d ",mu);
	if(nvparal_dir==0) master_printf("(none)");
	master_printf("\n");
	
	//assign coords of all ranks
	for(int loclx=0;loclx<loc_vol;loclx++)
	  {
	    coords loclx_vrank_coord; //host temporarily the vrank coord of loclx
	    for(int mu=0;mu<NDIM;mu++) loclx_vrank_coord[mu]=loc_coord_of_loclx[loclx][mu]/vloc_size[mu];
	    vrank_of_loclx[loclx]=lx_of_coord(loclx_vrank_coord,nvranks_per_dir);
	  }
	
	//offset between lx sites relative to a single vector
	for(int vrank=0;vrank<nvranks;vrank++)
	  {
	    coords c;
	    for(int mu=0;mu<NDIM;mu++) c[mu]=vrank_coord[vrank][mu]*vloc_size[mu];
	    vrank_loclx_offset[vrank]=loclx_of_coord(c);
	  }
      }
    }
    
    //! destructor
    void destroy()
    {
      if(inited)
	{
	  inited=false;
	  nissa_free(vrank_of_loclx);
	}
    }
    
    vranks_grid_t() : inited(false) {}
    ~vranks_grid_t() {destroy();}
  };
  
  //! structure to hold ordering
  template <class base_type> struct vranks_ord_t
  {
    //! initialization flag
    int inited;
    //! reference virtual grid
    vranks_grid_t<base_type> *vg;
    //! virtual local element corresponding to vrank 0 equivalent of loclx
    int *vloc_of_loclx;
    //! loclx corresponding to virtual vrank 0 element
    int *loclx_of_vloc;
    //! virtual local element (e/o) corresponding to vrank 0 equivalent of loclx
    int *veos_of_loclx;
    //! loclx corresponding to vrank 0 element
    int *loclx_of_veos[2];
    //! virtual local element (e/o) corresponding to vrank 0 equivalent of loceo
    int *veos_of_loceo[2];
    //! loclx corresponding to vrank 0 element
    int *loceo_of_veos[2];
    
    //! return the vn=0 equivalent
    int vn0lx_of_loclx(int loclx)
    {return loclx_of_vloc[vloc_of_loclx[loclx]];}
    
    //! remap from vir to lx
    template <class T,class VT> void virsome_remap_to_lx(T *out,VT *in)
    {
      static_assert(sizeof(T)*vg->nvranks==sizeof(VT),"vector type is not nvranks time the plain type");
      virsome_remap_to_lx_internal(out,in,sizeof(base_type),sizeof(T)/sizeof(base_type),vg->nvranks,vg->vrank_of_loclx,vg->vrank_loclx_offset,loclx_of_vloc);
    }
    
    //! remap from vir to lx
    template <class VT,class T> void lx_remap_to_virsome(VT *out,T *in)
    {
      static_assert(sizeof(T)*vg->nvranks==sizeof(VT),"vector type is not nvranks time the plain type");
      lx_remap_to_virsome_internal(out,in,sizeof(base_type),sizeof(T)/sizeof(base_type),vg->nvranks,vg->vrank_of_loclx,vloc_of_loclx);
    }
    
    //! init the geometry
    void init(vranks_grid_t<base_type> &ext_vg,std::function<void(vranks_ord_t*)>fill_vloc_of_loclx)
    {
      if(!inited)
      {
	//init and chain-init vg, copying its reference
  	inited=true;
	ext_vg.init();
	vg=&ext_vg;
	
	//allocate
	vloc_of_loclx=nissa_malloc("vloc_of_loclx",loc_vol,int);
	veos_of_loclx=nissa_malloc("veos_of_loclx",loc_vol,int);
	loclx_of_vloc=nissa_malloc("loclx_of_vloc",vg->vloc_vol,int);
	for(int par=0;par<2;par++)
	  {
	    loclx_of_veos[par]=nissa_malloc("loclx_of_veos",vg->vloc_volh,int);
	    loceo_of_veos[par]=nissa_malloc("loceo_of_veos",vg->vloc_volh,int);
	    veos_of_loceo[par]=nissa_malloc("veos_of_loceo",loc_volh,int);
	  }
	
	//fill all indices
	fill_vloc_of_loclx(this);
	
	//fill the rest
	for(int loclx=0;loclx<loc_vol;loclx++)
	  {
	    int par=loclx_parity[loclx];
	    int loceo=loceo_of_loclx[loclx];
	    int vloc=vloc_of_loclx[loclx];
	    int veos=vloc/2;
	    veos_of_loclx[loclx]=veos;
	    veos_of_loceo[par][loceo]=veos;
	    if(vg->vrank_of_loclx[loclx]==0)
	      {
		loclx_of_vloc[vloc]=loclx;
		loclx_of_veos[par][veos]=loclx;
		loceo_of_veos[par][veos]=loceo;
	      }
	  }
      }
    }
    
    //
    void destroy()
    {
      if(inited)
      {
  	inited=false;
	
   	nissa_free(vloc_of_loclx);
   	nissa_free(loclx_of_vloc);
	
   	nissa_free(veos_of_loclx);
	
  	for(int par=0;par<2;par++)
  	  {
	    nissa_free(loclx_of_veos[par]);
	    
  	    nissa_free(veos_of_loceo[par]);
  	    nissa_free(loceo_of_veos[par]);
	  }
      }
    }
    //! unset
    vranks_ord_t() : inited(false),vg(NULL) {}
    ~vranks_ord_t() {destroy();}
  };
  
  // void lx_conf_remap_to_virlx(vir_oct_su3 *out,quad_su3 *in);
  // void lx_conf_remap_to_virlx_blocked(vir_su3 *out,quad_su3 *in);
  // void virlx_conf_remap_to_lx(quad_su3 *out,vir_oct_su3 *in);
  // void lx_conf_remap_to_vireo(vir_oct_su3 **out,quad_su3 *in);
  // void lx_conf_remap_to_single_vireo(vir_single_oct_su3 **out,quad_su3 *in);
  // void eo_conf_remap_to_vireo(vir_oct_su3 **out,quad_su3 **in);
  // void eo_conf_remap_to_single_vireo(vir_single_oct_su3 **out,quad_su3 **in);
  
  // void lx_quad_su3_remap_to_virlx(vir_quad_su3 *out,quad_su3 *in);
  // void virlx_quad_su3_remap_to_lx(quad_su3 *out,vir_quad_su3 *in);
  // inline void lx_clover_term_t_remap_to_virlx(vir_clover_term_t *out,clover_term_t *in)
  // {lx_quad_su3_remap_to_virlx(out,in);}
  // inline void virlx_clover_term_t_remap_to_lx(clover_term_t *out,vir_clover_term_t *in)
  // {virlx_quad_su3_remap_to_lx(out,in);}
  
  // void evn_or_odd_quad_su3_remap_to_virevn_or_odd(vir_quad_su3 *out,quad_su3 *in,int par);
  // void virevn_or_odd_quad_su3_remap_to_evn_or_odd(quad_su3 *out,vir_quad_su3 *in,int par);
  // inline void evn_or_odd_clover_term_t_remap_to_virevn_or_odd(vir_clover_term_t *out,clover_term_t *in,int par)
  // {evn_or_odd_quad_su3_remap_to_virevn_or_odd(out,in,par);}
  // inline void virevn_or_odd_clover_term_t_remap_to_evn_or_odd(clover_term_t *out,vir_clover_term_t *in,int par)
  // {virevn_or_odd_quad_su3_remap_to_evn_or_odd(out,in,par);}
  
  // void virlx_spincolor_remap_to_lx(spincolor *out,vir_spincolor *in);
  // void virevn_or_odd_spincolor_remap_to_evn_or_odd(spincolor *out,vir_spincolor *in,int par);
  // void virlx_spincolor_128_remap_to_lx(spincolor_128 *out,vir_spincolor_128 *in);
  // void lx_spincolor_remap_to_virlx(vir_spincolor *out,spincolor *in);
  // void lx_spincolor_128_remap_to_virlx(vir_spincolor_128 *out,spincolor_128 *in);
  // void vireo_spincolor_remap_to_lx(spincolor *out,vir_spincolor **in);
  // void evn_or_odd_spincolor_remap_to_virevn_or_odd(vir_spincolor *out,spincolor *in,int par);
  // void virlx_clover_t_term_remap_to_lx(vir_clover_term_t *out,clover_term_t *in);
  
  // void vireo_color_remap_to_lx(color *out,vir_color **in);
  // void virevn_or_odd_color_remap_to_evn_or_odd(color *out,vir_color *in,int par);
  // void virevn_or_odd_single_color_remap_to_evn_or_odd(color *out,vir_single_color *in,int par);
  // void lx_spincolor_remap_to_vireo(vir_spincolor **out,spincolor *in);
  // void lx_color_remap_to_vireo(vir_color **out,color *in);
  // void lx_color_remap_to_single_vireo(vir_single_color **out,color *in);
  // void evn_or_odd_color_remap_to_virevn_or_odd(vir_color *out,color *in,int par);
  // void evn_or_odd_color_remap_to_single_virevn_or_odd(vir_single_color *out,color *in,int par);
  
  // void evn_or_odd_complex_vect_remap_to_virevn_or_odd(vir_complex *out,complex *in,int par,int vl);
  // inline void evn_or_odd_inv_clover_term_t_remap_to_virevn_or_odd(vir_inv_clover_term_t *out,inv_clover_term_t *in,int par)
  // {evn_or_odd_complex_vect_remap_to_virevn_or_odd((vir_complex*)out,(complex*)in,par,sizeof(inv_clover_term_t)/sizeof(complex));}
  
  void set_vranks_geometry();
  void unset_vranks_geometry();
  
  EXTERN_GEOMETRY_VIR vranks_grid_t<double> vdouble_grid;
  EXTERN_GEOMETRY_VIR vranks_grid_t<float> vfloat_grid;
  
  EXTERN_GEOMETRY_VIR vranks_ord_t<double> vlx_double_geom,vsf_double_geom;
  EXTERN_GEOMETRY_VIR vranks_ord_t<float> vlx_float_geom,vsf_float_geom;
  
#define DEFINE_VTYPES(VSHORT,LONG)					\
  typedef vranks_grid_t<LONG>::vcomplex NAME2(VSHORT,complex);		\
  typedef vranks_grid_t<LONG>::vcolor NAME2(VSHORT,color);		\
  typedef vranks_grid_t<LONG>::vsu3 NAME2(VSHORT,su3);			\
  typedef vranks_grid_t<LONG>::vhalfspincolor NAME2(VSHORT,halfspincolor); \
  typedef vranks_grid_t<LONG>::vcolor_halfspincolor NAME2(VSHORT,color_halfspincolor); \
  typedef vranks_grid_t<LONG>::vhalfspincolor_halfspincolor NAME2(VSHORT,halfspincolor_halfspincolor); \
  typedef vranks_grid_t<LONG>::vquad_su3 NAME2(VSHORT,quad_su3);	\
  typedef vranks_grid_t<LONG>::voct_su3 NAME2(VSHORT,oct_su3);		\
  typedef vranks_grid_t<LONG>::vspincolor NAME2(VSHORT,spincolor);	\
  typedef vranks_grid_t<LONG>::vhalfspin NAME2(VSHORT,halfspin);	\
  typedef vranks_grid_t<LONG>::vclover_term_t NAME2(VSHORT,clover_term_t); \
  typedef vranks_grid_t<LONG>::vinv_clover_term_t NAME2(VSHORT,inv_clover_term_t); \
  
  DEFINE_VTYPES(vd,double);
  DEFINE_VTYPES(vf,float);
  #undef DEFINE_VTYPES
  
  EXTERN_GEOMETRY_VIR int &vd_loc_vol INIT_TO(vdouble_grid.vloc_vol);
  EXTERN_GEOMETRY_VIR int &vf_loc_vol INIT_TO(vfloat_grid.vloc_vol);
}

#undef EXTERN_GEOMETRY_VIR
#undef INIT_TO

#endif
