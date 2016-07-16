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
#include "base/thread_macros.hpp"
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
  
  //! remap to and from lx
  void something_remap_to_virsome_internal(void *out,void *in,int vol,int size_per_site,int nel_per_site,int nvranks,int *vrank_of_locsite,int *idx);
  void virsome_remap_to_something_internal(void *out,void *in,int vol,int size_per_site,int nel_per_site,int nvranks,int *vrank_of_locsite,int *vrank_locsite_offset,int *idx);
  
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
    //! offset between corresponding lx elements in vrank 0
    int vrank_loclx_offset[nvranks];
    //! offset between corresponding eo elements in vrank 0
    int vrank_loceo_offset[nvranks];
    //! virtual rank hosting an lx local lattice element
    int *vrank_of_loclx;
    //! virtual rank hosting an eo local lattice element
    int *vrank_of_loceo[2];
    
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
	master_printf("Initializing vgrid for data of size %d bits, nvranks=%d\n",8*sizeof(base_type),nvranks);
	vrank_of_loclx=nissa_malloc("vrank_of_loclx",loc_vol,int);
	for(int par=0;par<2;par++) vrank_of_loceo[par]=nissa_malloc("vrank_of_loceo",loc_volh,int);
	
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
	
	bool something_found=false;
	int nparal_dir_min=RAND_MAX;
	int surf_min=RAND_MAX;
	while(vranks_partitioning.find_next_valid_partition(VRPD,loc_size,min_loc_size,fix_nvranks))
	  {
	    //set the temporary vir local size and compute number of parallel dirs
	    coords VLS;
	    int NP=0,SU=0;
	    for(int mu=0;mu<NDIM;mu++)
	      {
		VLS[mu]=loc_size[mu]/VRPD[mu];
		if(VRPD[mu]>1)
		  {
		    NP++;
		    SU+=vloc_vol/VLS[mu];
		  }
	      }
	    
	    //criterion: minimal number of directions, then minimal surf
	    bool criterion=false;
	    if(NP<nparal_dir_min)
	      {
		criterion=true;
		master_printf("Beating previous partitioning due to smaller number of virtually parallelized dirs (%d)\n",NP);
	      }
	    else
	      if(NP==nparal_dir_min&&SU<surf_min)
		{
		  criterion=true;
		  master_printf("Beating previous partitioning due to smaller virtual surface (%d)\n",SU);
		}
	    
	    //if fullfilling criterion
	    if(criterion)
	      {
		//store number of parallel dirs and surg
		nparal_dir_min=NP;
		surf_min=SU;
		//store number of ranks per dir
		for(int mu=0;mu<NDIM;mu++)
		  {
		    nvranks_per_dir[mu]=VRPD[mu];
		    vloc_size[mu]=loc_size[mu]/VRPD[mu];
		  }
		
		//print some infos
		master_printf(" new virtual local volume\t%d",vloc_size[0]);
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
	    vrank_of_loclx[loclx]=
	      vrank_of_loceo[loclx_parity[loclx]][loceo_of_loclx[loclx]]=
	      lx_of_coord(loclx_vrank_coord,nvranks_per_dir);
	  }
	
	//offset between lx sites relative to a single vector
	for(int vrank=0;vrank<nvranks;vrank++)
	  {
	    coords c;
	    for(int mu=0;mu<NDIM;mu++) c[mu]=vrank_coord[vrank][mu]*vloc_size[mu];
	    vrank_loclx_offset[vrank]=loclx_of_coord(c);
	    vrank_loceo_offset[vrank]=vrank_loclx_offset[vrank]/2;
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
	  for(int par=0;par<2;par++) nissa_free(vrank_of_loceo[par]);
	}
    }
    
    vranks_grid_t() : inited(false) {}
    ~vranks_grid_t() {destroy();}
  };
  
  //! flatten a vector type: T[2][2][nvranks] -> T[4][nvranks]
  template <typename T> struct flattened_vec_type
  {
    typedef BASE_TYPE(T) base_type;
    typedef base_type type[NBASE_EL(T)/vranks_grid_t<base_type>::nvranks][vranks_grid_t<base_type>::nvranks];
  };
#define FLATTENED_VEC_TYPE(T) typename flattened_vec_type<T>::type
  
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
    
    //! crash if not matching
    template <class T,class VT> void check_matching_size()
    {static_assert(sizeof(T)*vg->nvranks==sizeof(VT),"vector type is not nvranks time the plain type");}
    
    //////
    
    //! remap from vir to lx
    template <class T,class VT> void virsome_remap_to_lx(T *out,VT *in)
    {
      check_matching_size<T,VT>();
      virsome_remap_to_something_internal(out,in,loc_vol,sizeof(base_type),sizeof(T)/sizeof(base_type),vg->nvranks,vg->vrank_of_loclx,vg->vrank_loclx_offset,loclx_of_vloc);
    }
    
    //remap an lx vector to vir[some] layout
    template <class VT,class T> THREADABLE_FUNCTION_5ARG(something_remap_to_virsome, VT*,out, T*,in, int,vol, int*,vrank_of_locsite, int*,idx_out)
    {
      static_assert(NBASE_EL(VT)==NBASE_EL(T)*vg->nvranks,"number of vrank el does not match nvranks times the number of el");
      
      GET_THREAD_ID();
      //START_TIMING(remap_time,nremap);
      
      if((T*)out==in) CRASH("cannot use with out==in");
      master_printf("nbase_el: %d, %s\n",NBASE_EL(T),typeid(FLATTENED_TYPE(T)).name());
      master_printf("nbase_el v: %d, %s\n",NBASE_EL(VT),typeid(FLATTENED_TYPE(VT)).name());
      
      //copy the various virtual ranks
      NISSA_PARALLEL_LOOP(isite,0,vol)
	for(int iel=0;iel<NBASE_EL(T);iel++)
	  ((FLATTENED_VEC_TYPE(VT)*)out)[idx_out[isite]][vrank_of_locsite[isite]][iel]=
	    ((FLATTENED_TYPE(T)*)in)[isite][iel];
	   
      //wait filling
      set_borders_invalid(out);
      
      //STOP_TIMING(remap_time);
    }
    THREADABLE_FUNCTION_END
    
    //! remap from lx to vir
    template <class VT,class T> void lx_remap_to_virsome(VT *out,T *in)
    {
      //check_matching_size<T,VT>();
      something_remap_to_virsome(out,in,loc_vol,vg->vrank_of_loclx,vloc_of_loclx);
    }
    
    ///////
    
    //! remap from vir to evn or odd
    template <class T,class VT> void virsome_remap_to_evn_or_odd(T *out,VT *in,int par)
    {
      check_matching_size<T,VT>();
      virsome_remap_to_something_internal(out,in,loc_volh,sizeof(base_type),sizeof(T)/sizeof(base_type),vg->nvranks,vg->vrank_of_loceo[par],vg->vrank_loceo_offset,loceo_of_veos[par]);
    }
    
    //! remap from evn or odd to vir
    template <class VT,class T> void evn_or_odd_remap_to_virsome(VT *out,T *in,int par)
    {
      check_matching_size<T,VT>();
      something_remap_to_virsome_internal(out,in,loc_volh,sizeof(base_type),sizeof(T)/sizeof(base_type),vg->nvranks,vg->vrank_of_loceo[par],veos_of_loceo[par]);
    }
    
    ////////
    
    //! remap from vir to eo
    template <class T,class VT> void virsome_remap_to_eo(T *out[2],VT *in[2]) {for(int par=0;par<2;par++) virsome_remap_to_evn_or_odd(out[par],in[par],par);}
    
    //! remap from eo to vir
    template <class VT,class T> void eo_remap_to_virsome(VT *out[2],T *in[2]) {for(int par=0;par<2;par++) evn_or_odd_remap_to_virsome(out[par],in[par],par);}
    
    ////////
    
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
  
  /////////////////////////////////// automatic vectorization of type ///////////////////
  
  //! return false unless is a vrank type
  template <typename T> struct is_vrank_type
  {const static bool flag=false;};
#define IS_VRANK_TYPE(T) is_vrank_type<T>::flag
  
  //! specialize for really vectorized type (not automatic, unfortunately)
#define MARK_IS_VECTORIZED(TP)			\
  template<> struct is_vrank_type<TP>		\
  {const static bool flag=true;}
  
  //! vectorize a simple type
  template <typename _Tp> struct vectorize_type
  {typedef _Tp type[vranks_grid_t<_Tp>::nvranks];};
#define VECTORIZED_TYPE(TP) vectorize_type<TP>::type
  
  //! chain-vectorize arrays
  template <typename _Tp,std::size_t _Size> struct vectorize_type<_Tp[_Size]>
  {typedef typename VECTORIZED_TYPE(_Tp) type[_Size];};
  
  //! automatize the definition and marking
#define DEFINE_VECTORIZED_TYPE_NOMARK(TP)	\
  typedef VECTORIZED_TYPE(TP) NAME2(vd,TP)
#define DEFINE_VECTORIZED_TYPE(TP)			\
  DEFINE_VECTORIZED_TYPE_NOMARK(TP);			\
  MARK_IS_VECTORIZED(NAME2(vd,TP))
  
  DEFINE_VECTORIZED_TYPE(complex);
  DEFINE_VECTORIZED_TYPE(color);
  DEFINE_VECTORIZED_TYPE(su3);
  DEFINE_VECTORIZED_TYPE(halfspincolor);
  DEFINE_VECTORIZED_TYPE(color_halfspincolor);
  DEFINE_VECTORIZED_TYPE(halfspincolor_halfspincolor);
  DEFINE_VECTORIZED_TYPE(quad_su3);
  DEFINE_VECTORIZED_TYPE(oct_su3);
  DEFINE_VECTORIZED_TYPE(spincolor);
  DEFINE_VECTORIZED_TYPE(halfspin);
  DEFINE_VECTORIZED_TYPE_NOMARK(clover_term_t); //matches quad_su3
  DEFINE_VECTORIZED_TYPE(inv_clover_term_t);
  
  //////////////////// single version ////////////////////
  
  //! convert to single "any" (e.g. double) scalar type
  template <typename _Tp> struct float_type
  {typedef float type;};
#define FLOATED_TYPE(TP) float_type<TP>::type
  
  //! chain-vectorize conversion of arrays to float
  template <typename _Tp,std::size_t _Size> struct float_type<_Tp[_Size]>
  {typedef typename FLOATED_TYPE(_Tp) type[_Size];};
#define VECTORIZED_FLOATED_TYPE(TP) VECTORIZED_TYPE(FLOATED_TYPE(TP))
#define DEFINE_VECTORIZED_FLOATED_TYPE_NOMARK(TP)		\
  typedef VECTORIZED_FLOATED_TYPE(TP) NAME2(vf,TP)
#define DEFINE_VECTORIZED_FLOATED_TYPE(TP)			\
  DEFINE_VECTORIZED_FLOATED_TYPE_NOMARK(TP);			\
  MARK_IS_VECTORIZED(NAME2(vf,TP))
  
  DEFINE_VECTORIZED_FLOATED_TYPE(complex);
  DEFINE_VECTORIZED_FLOATED_TYPE(color);
  DEFINE_VECTORIZED_FLOATED_TYPE(su3);
  DEFINE_VECTORIZED_FLOATED_TYPE(halfspincolor);
  DEFINE_VECTORIZED_FLOATED_TYPE(color_halfspincolor);
  DEFINE_VECTORIZED_FLOATED_TYPE(halfspincolor_halfspincolor);
  DEFINE_VECTORIZED_FLOATED_TYPE(quad_su3);
  DEFINE_VECTORIZED_FLOATED_TYPE(oct_su3);
  DEFINE_VECTORIZED_FLOATED_TYPE(spincolor);
  DEFINE_VECTORIZED_FLOATED_TYPE(halfspin);
  DEFINE_VECTORIZED_FLOATED_TYPE_NOMARK(clover_term_t);
  DEFINE_VECTORIZED_FLOATED_TYPE(inv_clover_term_t);
  
  //////////////////// grids ///////////////////
  
  EXTERN_GEOMETRY_VIR int &vd_loc_vol INIT_TO(vdouble_grid.vloc_vol);
  EXTERN_GEOMETRY_VIR int &vf_loc_vol INIT_TO(vfloat_grid.vloc_vol);
  EXTERN_GEOMETRY_VIR int &vd_loc_volh INIT_TO(vdouble_grid.vloc_volh);
  EXTERN_GEOMETRY_VIR int &vf_loc_volh INIT_TO(vfloat_grid.vloc_volh);
}

#undef EXTERN_GEOMETRY_VIR
#undef INIT_TO

#endif
