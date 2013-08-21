#include "../../base/global_variables.h"
#include "../../base/thread_macros.h"
#include "../../new_types/new_types_definitions.h"
#include "../../routines/thread.h"

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
	int sour_gamma_comp=base_gamma[sour_igamma].entr[sour_gamma_sink_comp_id][RE]+
	  base_gamma[sour_igamma].entr[sour_gamma_sink_comp_id][IM];
	int sink_gamma_comp=base_gamma[sink_igamma].entr[sink_gamma_sink_comp_id][RE]+
	  base_gamma[sink_igamma].entr[sink_gamma_sink_comp_id][IM];
	
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
}

//print the content
void two_pts_comp_t::print(FILE *fout=stdout)
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
#ifndef USE_THREAD
  int start_contr_t=0,end_contr_t=loc_size[0]*ncontrib;
#else
  NISSA_CHUNK_WORKLOAD(start_contr_t,chunk_load_contr_t,end_contr_t,0,loc_size[0]*ncontrib,thread_id,NACTIVE_THREADS);
  chunk_load_contr_t++;
#endif
  
  //loop on the whole list
  int icontrib_t=0;
  double temp[ncontrib*loc_size[0]];
  for(int t=0;t<loc_size[0];t++)
    {
      double *S_forw_t=S_forw+t*nel*32;
      double *S_back_t=S_back+t*nel*32;
      
      //loop on all contractions
      for(two_pts_comp_t::iterator it=this->begin();it!=this->end();it++)
	{
	  if(icontrib_t>=start_contr_t&&icontrib_t<end_contr_t)
	    {
	      forw_back_comp_id_t comp_id(it->first);
	      corr_id_weight_t corr_id_weight_list(it->second);
	      
	      //compute the contribution and summ it to all the correlation where it contribute
	      double *S_forw_t_id=S_forw_t+nel*comp_id.first;
	      double *S_back_t_id=S_back_t+nel*comp_id.second;
	      
	      //compute the local product
	      double loc_temp=0;
	      for(int iel=0;iel<nel;iel++) loc_temp+=S_forw_t_id[iel]*S_back_t_id[iel];
	      temp[icontrib_t]=loc_temp;
	    }
	  icontrib_t++;
	}
    }

  THREAD_BARRIER();
  
  //summ the result to all the corr to which it contributes
  NISSA_PARALLEL_LOOP(loc_t,0,loc_size[0])
    {
      //set in and out
      icontrib_t=ncontrib*loc_size[0];
      int t=(glb_size[0]+loc_size[0]*rank_coord[0]+loc_t-twall)%glb_size[0];

      //loop on all contributions
      for(two_pts_comp_t::iterator it=this->begin();it!=this->end();it++)
	{
	  //loop on all correlation to which it contributes
	  for(corr_id_weight_t::iterator jt=it->second.begin();jt!=it->second.end();jt++)
	    {
	      uint16_t corr_id=jt->first;
	      double weight=jt->second;
	      
	      out[glb_size[0]*corr_id+t]+=weight*temp[icontrib_t];
	    }
	  icontrib_t++;
	}
    }
  
  THREAD_BARRIER();
}

//print optimized contractions to file
void print_optimized_contractions_to_file(FILE *fout,int ncontr,int *op_sour,int *op_sink,complex *contr,int twall,
					  const char *tag,double norm)
{
  if(norm!=1) crash("this would not be optimized");
  if(rank==0)
    for(int icontr=0;icontr<ncontr;icontr++)
      {
        fprintf(fout,"\n");
        fprintf(fout," # %s%s%s\n",tag,gtag[op_sink[icontr]],gtag[op_sour[icontr]]);
        for(int t=0;t<glb_size[0];t++) fprintf(fout,"%+016.16g\t%+016.16g\n",
					       ((double*)contr)[glb_size[0]*(0+2*icontr)+t]*norm,((double*)contr)[glb_size[0]*(1+2*icontr)+t]*norm);
      }
}
