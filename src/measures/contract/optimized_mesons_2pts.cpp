#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/dirac.hpp"
#include "threads/threads.hpp"

#include "optimized_mesons_2pts.hpp"

#ifdef BGQ
 #include "bgq/intrinsic.hpp"
#endif

namespace nissa
{
  //initialize the 2pts dirac matrix structure
  void two_pts_comp_t::add_sink_source_corr(uint16_t corr_id,double weight,int re_im,uint8_t sink_igamma,uint8_t sour_igamma)
  {
    //take pattern, according to the number of "i"
    int sink_pattern=(int)(base_gamma[sink_igamma].entr[0][IM]!=0);
    int sour_pattern=(int)(base_gamma[sour_igamma].entr[0][IM]!=0);
    int pattern=re_im+sink_pattern+sour_pattern; //count the number of "i"
    pattern_list.push_back(sink_pattern+2*(sour_pattern+2*re_im));
    for(int sour_gamma_sink_comp_id=0;sour_gamma_sink_comp_id<4;sour_gamma_sink_comp_id++)
      for(int sink_gamma_sink_comp_id=0;sink_gamma_sink_comp_id<4;sink_gamma_sink_comp_id++)
	{
	  //get components source id
	  int sour_gamma_sour_comp_id=base_gamma[sour_igamma].pos[sour_gamma_sink_comp_id];
	  int sink_gamma_sour_comp_id=base_gamma[sink_igamma].pos[sink_gamma_sink_comp_id];
	  
	  //get components value
	  int sour_gamma_comp=(int)(base_gamma[sour_igamma].entr[sour_gamma_sink_comp_id][RE]+
				    base_gamma[sour_igamma].entr[sour_gamma_sink_comp_id][IM]);
	  int sink_gamma_comp=int(base_gamma[sink_igamma].entr[sink_gamma_sink_comp_id][RE]+
				  base_gamma[sink_igamma].entr[sink_gamma_sink_comp_id][IM]);
	  
	  for(int re_im_S_forw=0;re_im_S_forw<2;re_im_S_forw++)
	    {
	      int re_im_S_pattern[2][2]={{0,1},{1,0}},re_im_S_back=re_im_S_pattern[pattern%2][re_im_S_forw];
	      
	      //find the entry of S to catch for each set
	      int forw_comp_id=re_im_S_forw+2*(sour_gamma_sink_comp_id+4*sink_gamma_sour_comp_id);
	      int back_comp_id=re_im_S_back+2*(sour_gamma_sour_comp_id+4*sink_gamma_sink_comp_id);
	      
	      //taking imaginary part is as taking Im(c)=-Re(i*c), so we count an "i" and put a -1
	      int re_im_sign[2]={+1,-1};
	      //gamma5 contains a - on 2 and 3 id of sour_gamma_sour_comp_id and sink_gamma_sink_comp_id
	      int g5_sign[4]={+1,+1,-1,-1};
	      //first contribution comes from the real part of forw prop
	      //first case is when doing re*re+im*im, second is re*im-im*re
	      int re_im_prod_sign[2][2]={{+1,+1},{+1,-1}};
	      //only if the two or three gammas are imaginary we have a -
	      int sour_sink_gamma_sign_list[4]={+1,+1,-1,-1},sour_sink_gamma_sign=sour_sink_gamma_sign_list[pattern];
	      
	      //compute the local sign
	      int loc_sign=re_im_prod_sign[pattern%2][re_im_S_forw]*re_im_sign[re_im]*
		g5_sign[sour_gamma_sour_comp_id]*g5_sign[sink_gamma_sink_comp_id]*
		sour_gamma_comp*sink_gamma_comp*sour_sink_gamma_sign;
	      
	      //add contribution to contraction
	      ((*this)[forw_back_comp_id_t(forw_comp_id,back_comp_id)])[corr_id]+=weight*loc_sign;
	    }
	}
    
    //scream zero elements
    this->scream();
  }
  
  //remove entries with weight 0
  void two_pts_comp_t::scream()
  {
    //erase ciclicly until nothing must be erased in it
    bool er_it;
    do
      {
	er_it=false;
	two_pts_comp_t::iterator it=this->begin();
	//loop on the whole list to search for something to be erased in it
	while(er_it==false&&it!=this->end())
	  {
	    //erase ciclicly until nothing must be erased in jt
	    bool er_jt;
	    do
	      {
		//mark no erasion and start from beginning
		er_jt=false;
		corr_id_weight_t::iterator jt=it->second.begin();
		
		//loop on the whole list to search for something to be erased
		while(er_jt==false&&jt!=it->second.end())
		  {
		    if(fabs(jt->second)<1.e-20) er_jt=true;
		    
		    //if not erasing, increase jt
		    if(!er_jt) jt++;
		  }
		
		//if erasing, erase jt
		if(er_jt==true) it->second.erase(jt);
	      }
	    while(er_jt);
	    
	    //check if to erease it
	    if(it->second.size()==0) er_it=true;
	    
	    //if not ereasing, increase it
	    if(!er_it) it++;
	  }
	
	//if erasing, erase it
	if(er_it) this->erase(it);
      }
    while(er_it);
  }
  
  //print the content
  void two_pts_comp_t::print(FILE *fout)
  {
    //loop on the whole list
    for(two_pts_comp_t::iterator it=this->begin();it!=this->end();it++)
      {
	forw_back_comp_id_t comp_id(it->first);
	corr_id_weight_t corr_id_weight_list(it->second);
	
	//intestation
	fprintf(fout,"%u %u: ",comp_id.forw(),comp_id.back());
	for(corr_id_weight_t::iterator jt=corr_id_weight_list.begin();jt!=corr_id_weight_list.end();jt++)
	  {
	    uint16_t corr_id=jt->first;
	    double weight=jt->second;
	    
	    fprintf(fout,"(%u,%+g) ",corr_id,weight);
	  }
	fprintf(fout,"\n");
      }
  }
  
  //summ the result (no glb reduction)
  void two_pts_comp_t::summ_the_loc_forw_back_contractions(double *out,double *S_forw,double *S_back,int nel,int twall)
  {
    GET_THREAD_ID();
    
    //define the workload (to be improved)
    int ncontrib=this->size();
#if THREADS_TYPE != OPENMP_THREADS
    int start_contr_t=0,end_contr_t=loc_size[0]*ncontrib;
#else
    NISSA_CHUNK_WORKLOAD(start_contr_t,chunk_load_contr_t,end_contr_t,0,loc_size[0]*ncontrib,thread_id,NACTIVE_THREADS);
    chunk_load_contr_t++;
#endif
    
    //loop on the whole list
    int icontrib_t=0;
    double *temp=nissa_malloc("temp",ncontrib*loc_size[0],double);
    for(int t=0;t<loc_size[0];t++)
      {
	double *S_forw_t=S_forw+t*nel*32;
	double *S_back_t=S_back+t*nel*32;
	
	//loop on all contractions
	for(two_pts_comp_t::iterator it=this->begin();it!=this->end();it++)
	  {
	    if(icontrib_t>=start_contr_t&&icontrib_t<end_contr_t)
	      {
		//compute the contribution and summ it to all the correlation where it contribute
		double *S_forw_t_id=S_forw_t+nel*it->first.first;
		double *S_back_t_id=S_back_t+nel*it->first.second;
		
		//compute the local product
#if defined BGQ && !defined BGQ_EMU
		S_forw_t_id-=4;
		S_back_t_id-=4;
		DECLARE_REG_VIR_COMPLEX(reg_loc_temp);
		REG_SPLAT_VIR_COMPLEX(reg_loc_temp,0);
		for(int iel=0;iel<nel;iel+=4)
		  {
		    DECLARE_REG_VIR_COMPLEX(reg_forw);
		    DECLARE_REG_VIR_COMPLEX(reg_back);
		    REG_LOAD_VIR_COMPLEX_AFTER_ADVANCING(reg_forw,S_forw_t_id);
		    REG_LOAD_VIR_COMPLEX_AFTER_ADVANCING(reg_back,S_back_t_id);
		    REG_VIR_COMPLEX_SUMM_THE_PROD_4DOUBLE(reg_loc_temp,reg_loc_temp,reg_forw,reg_back);
		  }
		
		double loc_temp[4];
		STORE_REG_VIR_COMPLEX(loc_temp,reg_loc_temp);
		temp[icontrib_t]=loc_temp[0]+loc_temp[1]+loc_temp[2]+loc_temp[3];
#else
		double loc_temp=0;
		for(int iel=0;iel<nel;iel++) loc_temp+=S_forw_t_id[iel]*S_back_t_id[iel];
		temp[icontrib_t]=loc_temp;
#endif
	      }
	    icontrib_t++;
	  }
      }
    
    THREAD_BARRIER();
    
    //summ the result to all the corr to which it contributes
    NISSA_PARALLEL_LOOP(loc_t,0,loc_size[0])
      {
	#warning not implemented
	// //set in and out
	// icontrib_t=ncontrib*loc_t;
	// int t=(glb_size[0]+loc_size[0]*rank_coord[0]+loc_t-twall)%glb_size[0];
	
	// //loop on all contributions
	// for(two_pts_comp_t::iterator it=this->begin();it!=this->end();it++)
	//   {
	//     //loop on all correlation to which it contributes
	//     for(corr_id_weight_t::iterator jt=it->second.begin();jt!=it->second.end();jt++)
	//       {
	// 	uint16_t corr_id=jt->first;
	// 	double weight=jt->second;
		
	// 	out[glb_size[0]*corr_id+t]+=weight*temp[icontrib_t];
	//       }
	//     icontrib_t++;
	//   }
      }
    NISSA_PARALLEL_LOOP_END;
    
    nissa_free(temp);
  }
  
  //print optimized correlations to file
  void two_pts_comp_t::print_correlations_to_file(FILE *fout,double *corr)
  {
    if(rank==0)
      for(int icorr=0;icorr<ncorr;icorr++)
	{
	  fprintf(fout,"\n");
	  fprintf(fout," # %s\n",corr_name[icorr].c_str());
	  for(int t=0;t<glb_size[0];t++) fprintf(fout,"%+16.16lg\n",corr[glb_size[0]*icorr+t]);
	}
  }
  
  //summassign
  two_pts_comp_t operator+=(two_pts_comp_t &a,two_pts_comp_t &b)
  {
    for(two_pts_comp_t::iterator it=b.begin();it!=b.end();it++)
      for(corr_id_weight_t::iterator jt=it->second.begin();jt!=it->second.end();jt++)
	a[it->first][jt->first]+=jt->second;
    
    a.scream();
    
    return a;
  }
  
  //subtassign
  two_pts_comp_t operator-=(two_pts_comp_t &a,two_pts_comp_t &b)
  {
    for(two_pts_comp_t::iterator it=b.begin();it!=b.end();it++)
      for(corr_id_weight_t::iterator jt=it->second.begin();jt!=it->second.end();jt++)
	a[it->first][jt->first]-=jt->second;
    
    a.scream();
    
    return a;
  }
  
  //prodassign with double
  two_pts_comp_t operator*=(two_pts_comp_t &a,double b)
  {
    for(two_pts_comp_t::iterator it=a.begin();it!=a.end();it++)
      for(corr_id_weight_t::iterator jt=it->second.begin();jt!=it->second.end();jt++)
	(jt->second)*=b;
    
    a.scream();
    
    return a;
  }
  
  //divassign with double
  two_pts_comp_t operator/=(two_pts_comp_t &a,double b)
  {
    for(two_pts_comp_t::iterator it=a.begin();it!=a.end();it++)
      for(corr_id_weight_t::iterator jt=it->second.begin();jt!=it->second.end();jt++)
	(jt->second)/=b;
    
    a.scream();
    
    return a;
  }
}
