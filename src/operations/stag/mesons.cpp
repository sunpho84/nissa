#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"
#include "operations/gauge_fixing.hpp"

#include "stag.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "mesons.hpp"

//in this formalism, shift=t+2*(x+2*(y+2*z))

namespace nissa
{
  namespace
  {
    int nop;
    int ncombo;
    int nflavs;
    
    //form the mask for x (-1)^[x*(s^<+n^>)](-1)^[n*(s+n)^<]
    int form_stag_meson_pattern(int ispin,int itaste)
    {
      //add g5*g5
      ispin^=15;
      itaste^=15;
      
      int res=0;
      for(int mu=0;mu<NDIM;mu++)
	{
	  int p=0;
	  for(int nu=0;nu<mu;nu++) p+=(itaste>>nu)&1;
	  for(int nu=mu+1;nu<NDIM;nu++) p+=(ispin>>nu)&1;
	  p&=1;
	  
	  res+=(p<<mu);
	}
      
      return res;
    }
    
    //compute the index where to store
    int icombo(int iflav,int iop,int t)
    {return t+glb_size[0]*(iop+nop*iflav);}
    
    inline void addrem_stagphases(quad_su3 **conf)
    {
      GET_THREAD_ID();
      
      for(int eo=0;eo<2;eo++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      coords ph;
	      get_stagphase_of_lx(ph,loclx_of_loceo[eo][ieo]);
	      for(int mu=0;mu<NDIM;mu++) su3_prodassign_double(conf[eo][ieo][mu],ph[mu]);
	    }
	  set_borders_invalid(conf[eo]);
	}
    }
    
    //apply the operator
    inline void apply_op(color **source,color **temp,quad_su3 **conf,int shift,color **ori_source)
    {
      GET_THREAD_ID();
      
      for(int eo=0;eo<2;eo++) vector_copy(source[eo],ori_source[eo]);
      
      addrem_stagphases(conf);
      
      for(int mu=0;mu<NDIM;mu++)
	if((shift>>mu)&0x1)
	  {
	    //write header, copy and communicate
	    verbosity_lv2_master_printf("shift %d %d\n",shift,mu);
	    for(int eo=0;eo<2;eo++) vector_copy(temp[eo],source[eo]);
	    communicate_ev_and_od_color_borders(temp);
	    communicate_ev_and_od_quad_su3_borders(conf);
	    
	    for(int eo=0;eo<2;eo++)
	      NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
		{
		  int up=loceo_neighup[eo][ieo][mu];
		  int dw=loceo_neighdw[eo][ieo][mu];
		  unsafe_su3_prod_color(source[eo][ieo],conf[eo][ieo][mu],temp[!eo][up]);
		  su3_dag_summ_the_prod_color(source[eo][ieo],conf[!eo][dw][mu],temp[!eo][dw]);
		  color_prod_double(source[eo][ieo],source[eo][ieo],0.5);
		}
	    for(int eo=0;eo<2;eo++) set_borders_invalid(source[eo]);
	  }
      
      addrem_stagphases(conf);
    }
    
    //add the phases
    inline void put_phases(color **source,int mask)
    {
      GET_THREAD_ID();
      
      //put the phases
      for(int eo=0;eo<2;eo++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      int sign=1,ivol=loclx_of_loceo[eo][ieo];
	      for(int mu=0;mu<NDIM;mu++) sign*=1-2*(((mask>>mu)&0x1)&&(glb_coord_of_loclx[ivol][mu]&0x1));
	      if(abs(sign)!=1) crash("unexpected sign %d",sign);
	      color_prod_double(source[eo][ieo],source[eo][ieo],sign);
	    }
	  set_borders_invalid(source[eo]);
	}
    }
    
    //check it is an hypercube origin
    inline int is_hypercube_shift(int ivol,int ishift)
    {
      int is=1;
      for(int mu=0;mu<NDIM;mu++) is&=((glb_coord_of_loclx[ivol][mu]%2)==((ishift>>mu)%2));
      return is;
    }
  }
  
  //compute correlation functions for staggered mesons, arbitary taste and spin
  THREADABLE_FUNCTION_4ARG(compute_staggered_meson_corr, complex*,corr, quad_su3**,conf, theory_pars_t*,tp, stag_meson_corr_meas_pars_t*,meas_pars)
  {
    GET_THREAD_ID();
    
    //allocate
    color *ori_source[2],*source[2],*sol[2],*quark[nop][2];
    for(int eo=0;eo<2;eo++) ori_source[eo]=nissa_malloc("ori_source",loc_volh+bord_volh,color);
    for(int eo=0;eo<2;eo++) source[eo]=nissa_malloc("source",loc_volh+bord_volh,color);
    for(int eo=0;eo<2;eo++) sol[eo]=nissa_malloc("sol",loc_volh+bord_volh,color);
    for(int iop=0;iop<nop;iop++)
      for(int eo=0;eo<2;eo++)
	quark[iop][eo]=nissa_malloc("quark",loc_volh+bord_volh,color);
    complex *loc_corr=new complex[ncombo];
    memset(loc_corr,0,sizeof(complex)*ncombo);
    
    //form the masks
    int mask[nop],shift[nop];
    for(int iop=0;iop<nop;iop++)
      {
	int spin=meas_pars->mesons[iop].first;
	int taste=meas_pars->mesons[iop].second;
	shift[iop]=(spin^taste);
	mask[iop]=form_stag_meson_pattern(spin,taste);
	if((shift[iop])&1) crash("operator %d (%d %d) has unpaired number of g0",iop,spin,taste);
	master_printf(" iop %d (%d %d),\tmask: %d,\tshift: %d\n",iop,spin,taste,mask[iop],shift[iop]);
      }
    
    for(int ihit=0;ihit<meas_pars->nhits;ihit++)
      {
	generate_fully_undiluted_eo_source(ori_source,RND_Z4,0);
	for(int eo=0;eo<2;eo++)
	  {
	    NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	      if(!is_hypercube_shift(loclx_of_loceo[eo][ieo],0)) color_put_to_zero(ori_source[eo][ieo]);
	    set_borders_invalid(ori_source[eo]);
	  }
	
	for(int iflav=0;iflav<nflavs;iflav++)
	  {
	    for(int iop=0;iop<nop;iop++)
	      {
		apply_op(source,sol,conf,shift[iop],ori_source);
		put_phases(source,mask[iop]);
		mult_Minv(sol,conf,tp,iflav,meas_pars->residue,source);
		put_phases(sol,mask[iop]);
		apply_op(quark[iop],source,conf,shift[iop],sol);
	      }
	    
	    for(int iop=0;iop<nop;iop++)
	      for(int eo=0;eo<2;eo++)
		NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
		  {
		    int ivol=loclx_of_loceo[eo][ieo];
		    int t=glb_coord_of_loclx[ivol][0];
		    for(int ic=0;ic<NCOL;ic++)
		      complex_summ_the_conj2_prod(loc_corr[icombo(iflav,iop,t)],quark[0][eo][ieo][ic],quark[iop][eo][ieo][ic]);
		  }
	    THREAD_BARRIER();
	  }
      }
    
    //reduce
    glb_threads_reduce_double_vect((double*)loc_corr,2*ncombo);
    if(IS_MASTER_THREAD) glb_nodes_reduce_complex_vect(corr,loc_corr,ncombo);
    
    for(int eo=0;eo<2;eo++)
      {
	nissa_free(ori_source[eo]);
	nissa_free(source[eo]);
	nissa_free(sol[eo]);
      }
    for(int iop=0;iop<nop;iop++)
      for(int eo=0;eo<2;eo++)
	nissa_free(quark[iop][eo]);
    delete[] loc_corr;
  }
  THREADABLE_FUNCTION_END
  
  //compute and print
  void measure_staggered_meson_corr(quad_su3 **ext_conf,theory_pars_t &tp,stag_meson_corr_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    nop=meas_pars.mesons.size();
    nflavs=tp.nflavs;
    ncombo=icombo(nflavs-1,nop-1,glb_size[0]-1)+1;
    complex *corr=nissa_malloc("corr",ncombo,complex);
    
    compute_staggered_meson_corr(corr,ext_conf,&tp,&meas_pars);
    
    //open the file, allocate point result and source
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    for(int iop=0;iop<nop;iop++)
      for(int iflav=0;iflav<nflavs;iflav++)
	{
	  int spin=meas_pars.mesons[iop].first;
	  int taste=meas_pars.mesons[iop].second;
	  master_fprintf(file," # conf %d ; iop %d , spin %d , taste %d ; flv = %d , m = %lg\n",
			 iconf,iop,spin,taste,iflav,tp.quark_content[iflav].mass);
	  for(int t=0;t<glb_size[0];t++)
	    {
	      int ic=icombo(iflav,iop,t);
	      master_fprintf(file,"%d %+016.16lg %+016.016lg\n",t,corr[ic][RE],corr[ic][IM]);
	    }
	    master_fprintf(file,"\n");
	  }
    close_file(file);
    
    nissa_free(corr);
  }
  
  //nucleon correlators
  int stag_meson_corr_meas_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    if(flag||full)
      {
	nprinted+=nissa::master_fprintf(fout,"StagMesonCorrelators\n");
	if(flag!=1||full) nprinted+=nissa::master_fprintf(fout,"Each\t\t=\t%d\n",flag);
	if(path!=def_path()||full) nprinted+=nissa::master_fprintf(fout,"Path\t\t=\t\"%s\"\n",path.c_str());
	if(mesons.size()||full)
	  {
	    nprinted+=nissa::master_fprintf(fout,"Mesons\t\t=\t{");
	    for(size_t i=0;i<mesons.size();i++)
	      {
		nprinted+=nissa::master_fprintf(fout,"(%d,%d)",mesons[i].first,mesons[i].second);
		if(i!=mesons.size()-1) nprinted+=nissa::master_fprintf(fout,",");
		else                   nprinted+=nissa::master_fprintf(fout,"}\n");
	      }
	  }
	if(residue!=def_residue()||full) nprinted+=nissa::master_fprintf(fout,"Residue\t\t=\t\"%lg\"\n",residue);
	if(nhits!=def_nhits()||full) nprinted+=nissa::master_fprintf(fout,"NHits\t\t=\t%d\n",nhits);
      }
    else if(full) nprinted+=nissa::master_fprintf(fout,"StagMesonCorrelators No\n");
    
    return nprinted;
  }
}
