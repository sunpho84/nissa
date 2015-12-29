#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"
#include "operations/gauge_fixing.hpp"

#include "stag.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "mesons.hpp"

namespace nissa
{
  namespace
  {
    int nop;
    int ncombo;
    int nflavs;
    
    //form the mask for x (-1)^[x*(s^<+n^>)]
    inline int form_stag_meson_pattern(int ispin,int itaste)
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
    inline int icombo(int iflav,int iop,int t)
    {return t+glb_size[0]*(iop+nop*iflav);}
    
    //apply a single shift
    void apply_covariant_shift(color **out,quad_su3 **conf,int mu,color **in)
    {
      GET_THREAD_ID();
      
      communicate_ev_and_od_color_borders(in);
      communicate_ev_and_od_quad_su3_borders(conf);
      
      for(int eo=0;eo<2;eo++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      int up=loceo_neighup[eo][ieo][mu];
	      int dw=loceo_neighdw[eo][ieo][mu];
	      unsafe_su3_prod_color(out[eo][ieo],conf[eo][ieo][mu],in[!eo][up]);
	      su3_dag_summ_the_prod_color(out[eo][ieo],conf[!eo][dw][mu],in[!eo][dw]);
	      color_prod_double(out[eo][ieo],out[eo][ieo],0.5);
	    }
	  set_borders_invalid(out[eo]);
	}
    }
    
    //apply the operator
    inline void apply_op_single_perm(color **out,color **temp,quad_su3 **conf,std::vector<int> &list_dir,color **in)
    {
      //make a temporary copy
      for(int eo=0;eo<2;eo++) vector_copy(temp[eo],in[eo]);
      
      for(std::vector<int>::iterator mu_it=list_dir.begin();mu_it!=list_dir.end();mu_it++)
	{
	  //write header, copy and communicate
	  verbosity_lv2_master_printf(" shift %d\n",*mu_it);
	  if(mu_it!=list_dir.begin()) for(int eo=0;eo<2;eo++) vector_copy(temp[eo],out[eo]);
	  
	  //make the shift
	  apply_covariant_shift(out,conf,*mu_it,temp);
	}
    }
    
    //apply the operator summing all permutations
    inline void apply_op(color **out,color **single_perm,color **internal_temp,quad_su3 **conf,int shift,color **in)
    {
      //make a list that can be easily permuted
      std::vector<int> list_dir;
      for(int mu=0;mu<NDIM;mu++)
	if((shift>>mu)&0x1)
	  list_dir.push_back(mu);
      std::sort(list_dir.begin(),list_dir.end());
      
      if(list_dir.size())
	{
	  //summ all perms
	  int nperm=0;
	  for(int eo=0;eo<2;eo++) vector_reset(out[eo]);
	  do
	    {
	      //incrementing the number of permutations
	      verbosity_lv2_master_printf("Considering permutation %d:",nperm);
	      for(std::vector<int>::iterator it=list_dir.begin();it!=list_dir.end();it++) verbosity_lv2_master_printf(" %d",*it);
	      verbosity_lv2_master_printf("\n");
	      nperm++;
	      
	      //apply and summ
	      apply_op_single_perm(single_perm,internal_temp,conf,list_dir,in);
	      for(int eo=0;eo<2;eo++) double_vector_summassign((double*)(out[eo]),(double*)(single_perm[eo]),loc_volh*sizeof(color)/sizeof(double));
	    }
	  while(std::next_permutation(list_dir.begin(),list_dir.end()));
	  
	  //final normalization
	  for(int eo=0;eo<2;eo++) double_vector_prod_double((double*)(out[eo]),(double*)(out[eo]),1.0/nperm,loc_volh*sizeof(color)/sizeof(double));
	}
      else for(int eo=0;eo<2;eo++) vector_copy(out[eo],in[eo]);
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
    color *ori_source[2],*source[2],*sol[2],*quark[nop][2],*temp[2][2];
    for(int eo=0;eo<2;eo++) ori_source[eo]=nissa_malloc("ori_source",loc_volh+bord_volh,color);
    for(int eo=0;eo<2;eo++) source[eo]=nissa_malloc("source",loc_volh+bord_volh,color);
    for(int eo=0;eo<2;eo++) sol[eo]=nissa_malloc("sol",loc_volh+bord_volh,color);
    for(int iop=0;iop<nop;iop++)
      for(int eo=0;eo<2;eo++)
	quark[iop][eo]=nissa_malloc("quark",loc_volh+bord_volh,color);
    for(int itemp=0;itemp<2;itemp++)
      for(int eo=0;eo<2;eo++)
	temp[itemp][eo]=nissa_malloc("temp",loc_volh+bord_volh,color);
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
	if((shift[iop])&1) crash("operator %d (%d %d) has unmarched number of g0",iop,spin,taste);
	master_printf(" iop %d (%d %d),\tmask: %d,\tshift: %d\n",iop,spin,taste,mask[iop],shift[iop]);
      }
    
    for(int ihit=0;ihit<meas_pars->nhits;ihit++)
      {
	//generate tso
	int tso;
	if(IS_MASTER_THREAD) tso=rnd_get_unif(&glb_rnd_gen,0,glb_size[0]);
	THREAD_BROADCAST(tso,tso);
	master_printf("tsource: %d\n",tso);
	
	//generate source
	generate_fully_undiluted_eo_source(ori_source,RND_Z4,tso);
	
	for(int iflav=0;iflav<nflavs;iflav++)
	  {
	    for(int iop=0;iop<nop;iop++)
	      {
		apply_op(source,temp[0],temp[1],conf,shift[iop],ori_source);
		put_phases(source,mask[iop]);
		mult_Minv(sol,conf,tp,iflav,meas_pars->residue,source);
		apply_op(quark[iop],temp[0],temp[1],conf,shift[iop],sol);
		put_phases(quark[iop],mask[iop]);
	      }
	    
	    for(int iop=0;iop<nop;iop++)
	      for(int eo=0;eo<2;eo++)
		NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
		  {
		    int ivol=loclx_of_loceo[eo][ieo];
		    int t=(glb_coord_of_loclx[ivol][0]-tso+glb_size[0])%glb_size[0];
		    for(int ic=0;ic<NCOL;ic++)
		      complex_summ_the_conj1_prod(loc_corr[icombo(iflav,iop,t)],quark[0][eo][ieo][ic],quark[iop][eo][ieo][ic]);
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
    for(int itemp=0;itemp<2;itemp++)
      for(int eo=0;eo<2;eo++)
	nissa_free(temp[itemp][eo]);
    delete[] loc_corr;
  }
  THREADABLE_FUNCTION_END
  
  //compute and print
  void measure_staggered_meson_corr(quad_su3 **ext_conf,theory_pars_t &tp,stag_meson_corr_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    nop=meas_pars.mesons.size();
    nflavs=tp.nflavs;
    ncombo=icombo(nflavs-1,nop-1,glb_size[0]-1)+1;
    double norm=1.0/(meas_pars.nhits*glb_spat_vol);
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
	      master_fprintf(file,"%d %+016.16lg %+016.016lg\n",t,corr[ic][RE]*norm,corr[ic][IM]*norm);
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
