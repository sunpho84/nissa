#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "free_theory/free_theory_types.hpp"
#include "free_theory/tlSym_gauge_propagator.hpp"
#include "geometry/geometry_mix.hpp"
#include "measures/fermions/qed_corr.hpp"
#include "measures/fermions/stag.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

#include "mesons.hpp"

namespace nissa
{
  namespace stag
  {
    //return directly a eosplit photon field
    void get_eo_photon(eo_ptr<spin1field> out,gauge_info photon)
    {
      crash("reimplement");
      // //allocate lx version of photon field
      // spin1field *photon_eta=nissa_malloc("photon_eta",locVol+bord_vol,spin1field);
      // spin1field *photon_field=nissa_malloc("photon_field",locVol+bord_vol,spin1field);
      // spin1field *photon_phi=nissa_malloc("photon_phi",locVol+bord_vol,spin1field);
      
      // //generate source and stochastich propagator
      // generate_stochastic_tlSym_gauge_propagator(photon_phi,photon_eta,photon);
      // multiply_by_sqrt_tlSym_gauge_propagator(photon_field,photon_eta,photon);
      // split_lx_vector_into_eo_parts(out,photon_field);
      
      // nissa_free(photon_phi);
      // nissa_free(photon_field);
      // nissa_free(photon_eta);
    }
    
    void insert_tadpole_handle(complex out,eo_ptr<spin1field> aux,int par,int ieo,int mu,void *pars){out[RE]=(*(Momentum*)pars)[mu];out[IM]=0;}
    void insert_conserved_current_handle(complex out,eo_ptr<spin1field> aux,int par,int ieo,int mu,void *pars){out[RE]=((int*)pars)[mu];out[IM]=0;}
    void insert_time_conserved_vector_current_handle(complex out,eo_ptr<spin1field> aux,int par,int ieo,int mu,void *pars){out[RE]=(mu==0);out[IM]=0;}
    
    //insert the tadpol
    void insert_tadpole(eo_ptr<color> out,eo_ptr<quad_su3> conf,theory_pars_t* theory_pars,int iflav,eo_ptr<color> in,Momentum& tad,int t)
    {
      //call with no source insertion, plus between fw and bw, and a global -0.25
      complex fw_factor={-0.25,0},bw_factor={-0.25,0};
      insert_vector_vertex(out,conf,theory_pars,iflav,{NULL,NULL},in,fw_factor,bw_factor,insert_tadpole_handle,t,&tad);
    }
    
    //insert the external source, that is one of the two extrema of the stoch prop
    void insert_external_source(eo_ptr<color> out,eo_ptr<quad_su3> conf,theory_pars_t* theory_pars,int iflav,eo_ptr<spin1field> curr,eo_ptr<color> in,int t)
    {
      //call with source insertion, minus between fw and bw, and a global i*0.5
      complex fw_factor={0,+0.5},bw_factor={0,-0.5};
      insert_vector_vertex(out,conf,theory_pars,iflav,curr,in,fw_factor,bw_factor,insert_external_source_handle,t);
    }
    
    //insert the time componente of the vectorial current
    void insert_time_conserved_vector_current(eo_ptr<color> out,eo_ptr<quad_su3> conf,theory_pars_t* theory_pars,int iflav,eo_ptr<color> in,int t)
    {
      //call with no source insertion, minus between fw and bw, and a global i*0.5
      complex fw_factor={0,+0.5},bw_factor={0,-0.5};
      insert_vector_vertex(out,conf,theory_pars,iflav,{NULL,NULL},in,fw_factor,bw_factor,insert_time_conserved_vector_current_handle,t);
    }
  }
  
  using namespace stag;
  
  namespace
  {
    //const int nop_t=4;
    enum ins_t{S,T,F,V};
    //char op_name[nop_t][2]={"S","T","F","V"};
    //const int nprop_t=6;
    enum prop_t{P0,PS,PT,P1,P2,PV};
    //char prop_name[nprop_t][2]={"0","S","T","1","2","V"};
    
    struct ins_map_t
    {
      int sou;
      int op;
      ins_map_t(int sou,int op) : sou(sou),op(op) {}
      ins_map_t(){};
    };
    
  }
  
  //compute and print
  void measure_qed_corr(eo_ptr<quad_su3> conf,theory_pars_t theory_pars,qed_corr_meas_pars_t meas_pars,int iconf,int conf_created)
  {
    crash("reimplement");
    
    // //open the file, allocate point result and source
    // FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    // NEW_FIELD_T(ori_source);
    // NEW_FIELD_T(temp_source);
    
    // //set photon
    // gauge_info photon;
    // photon.alpha=FEYNMAN_ALPHA;
    // for(int mu=0;mu<NDIM;mu++) photon.bc[mu]=0;
    // photon.c1=C1_WILSON;
    // photon.zms=UNNO_ALEMANNA;
    crash(" ");
    
    // //compute tadpole
    // Momentum tadpole=compute_tadpole(photon);
    
    // //allocate
    // const int nflavs=theory_pars.nflavs();
    // eo_ptr<spin1field> photon_field={nissa_malloc("photon_phi_ev",locVolh+bord_volh,spin1field),nissa_malloc("photon_phi_od",locVolh+bord_volh,spin1field)};
    // eo_ptr<color> M[nprop_t*nflavs];
    // for(int i=0;i<nprop_t*nflavs;i++)
    //   for(int par=0;par<2;par++)
    // 	M[i][par]=nissa_malloc(combine("M_%d_%d",i,par).c_str(),locVolh+bord_volh,color);
    
    // //write the map of how to build props
    // std::vector<ins_map_t> prop_build(nprop_t);
    // prop_build[P0]=ins_map_t(-1,S);
    // prop_build[P1]=ins_map_t(P0,F);
    // prop_build[P2]=ins_map_t(P1,F);
    // prop_build[PS]=ins_map_t(P0,S);
    // prop_build[PT]=ins_map_t(P0,T);
    // prop_build[PV]=ins_map_t(P0,V);
    
    // //write how to meake the contractions
    // std::vector<std::pair<int,int> > contr_map;
    // contr_map.push_back(std::make_pair(P0,P0));
    // contr_map.push_back(std::make_pair(P0,PS));
    // contr_map.push_back(std::make_pair(P0,PT));
    // contr_map.push_back(std::make_pair(P1,P1));
    // contr_map.push_back(std::make_pair(P0,P2));
    // contr_map.push_back(std::make_pair(P0,PV));
    
    // //init the contr
    // //int ncontr_tot=contr_map.size()*nflavs*nflavs,contr_tot_size=ncontr_tot*glb_size[0];
    // complex *glb_contr=nullptr;
    // crash("#warning reimplement nissa_malloc(\"glb_contr\",contr_tot_size*nthreads,complex);");
    // // complex *loc_contr=glb_contr+THREAD_ID*contr_tot_size;
    
    // for(int icopy=0;icopy<meas_pars.ncopies;icopy++)
    //   {
    // 	vector_reset(glb_contr);
	
    // 	for(int ihit=0;ihit<meas_pars.nhits;ihit++)
    // 	  {
    // 	    verbosity_lv1_master_printf("Computing hit %d/%d\n",ihit,meas_pars.nhits);
	    
    // 	    //get global time
    // 	    int tso=rnd_get_unif(&glb_rnd_gen,0,glbSize[0]);
    // 	    verbosity_lv1_master_printf("tsource: %d\n",tso);
	    
    // 	    //generate sources
    // 	    get_eo_photon(photon_field,photon);
    // 	    fill_source(ori_source,tso,meas_pars.rnd_type);
	    
    // 	    // tso=0;
    // 	    // vector_reset(ori_source[EVN]);
    // 	    // vector_reset(ori_source[ODD]);
    // 	    // for(int icol=0;icol<3;icol++) ori_source[EVN][0][icol][RE]=1;
	    
    // 	    for(int iflav=0;iflav<nflavs;iflav++)
    // 	      for(size_t iprop=0;iprop<prop_build.size();iprop++)
    // 		{
    // 		  //select the source
    // 		  eo_ptr<color> so;
    // 		  if(prop_build[iprop].sou==-1) so=ori_source;
    // 		  else so=M[prop_build[iprop].sou+nprop_t*iflav];
		  
    // 		  //make the insertion
    // 		  verbosity_lv1_master_printf("Producing prop for flav %d, type %s, inserting operator %s on top of %s\n",
    // 					      iflav,prop_name[iprop],op_name[prop_build[iprop].op],
    // 					      (prop_build[iprop].sou==-1)?"so":prop_name[prop_build[iprop].sou]);
		  
    // 		  switch(prop_build[iprop].op)
    // 		    {
    // 		    case S:for(int par=0;par<2;par++) vector_copy(temp_source[par],so[par]);break;
    // 		    case T:insert_tadpole(temp_source,conf,&theory_pars,iflav,so,tadpole,-1);break;
    // 		    case F:insert_external_source(temp_source,conf,&theory_pars,iflav,photon_field,so,-1);break;
    // 		    case V:insert_time_conserved_vector_current(temp_source,conf,&theory_pars,iflav,so,(tso+glbSize[0]/4)%glbSize[0]);break;
    // 		    }
		  
    // 		  //invert
    // 		  // if(prop_build[iprop].op!=V)
    // 		    MINV(M[iprop+nprop_t*iflav],iflav,temp_source);
    // 		    //else for(int par=0;par<2;par++) vector_copy(M[iprop+nprop_t*iflav][par],temp_source[par]);
    // 		}
	    
    // 	    for(int iflav=0;iflav<nflavs;iflav++)
    // 	      for(int jflav=0;jflav<nflavs;jflav++)
    // 		for(size_t icontr=0;icontr<contr_map.size();icontr++)
    // 		  {
    // 		    // eo_ptr<color> A=(contr_map[icontr].first==-1)?ori_source:M[contr_map[icontr].first+nprop_t*iflav];
    // 		    // eo_ptr<color> B=(contr_map[icontr].second==-1)?ori_source:M[contr_map[icontr].second+nprop_t*jflav];
		    
    // 			  crash("#warning reimplement");
    // 		    for(int par=0;par<2;par++)
    // 		      NISSA_PARALLEL_LOOP(ieo,0,locVolh)
    // 			{
    // 			  // int ivol=loclx_of_loceo[par][ieo];
    // 			  // int t=(glb_coord_of_loclx[ivol][0]+glb_size[0]-tso)%glb_size[0];
    // 			  // for(int ic=0;ic<NCOL;ic++)
    // 			  //   complex_summ_the_conj1_prod(loc_contr[t+glb_size[0]*(icontr+contr_map.size()*(iflav+nflavs*jflav))],
    // 			  // 				A[par][ieo][ic],B[par][ieo][ic]);
    // 			}
    // 		    NISSA_PARALLEL_LOOP_END;
    // 		  }
    // 	  }
	
    // 	//reduce
    // 	crash("#warning reimplement glb_threads_reduce_complex_vect(loc_contr,contr_tot_size);");
    // 	crash("#warning if(IS_MASTER_THREAD) glb_nodes_reduce_complex_vect(glb_contr,contr_tot_size);");
	
    // 	//print
    // 	double norm=1.0/(meas_pars.nhits*glbSpatVol);
    // 	for(int iflav=0;iflav<nflavs;iflav++)
    // 	  for(int jflav=0;jflav<nflavs;jflav++)
    // 	    {
    // 	      master_fprintf(file,"\n # conf %d , iq_rev = %d , mq_rev = %lg , iq_ins = %d , mq_ins = %lg\n",
    // 			     iconf,iflav,theory_pars.quarks[iflav].mass,jflav,theory_pars.quarks[jflav].mass);
    // 	      for(size_t icontr=0;icontr<contr_map.size();icontr++)
    // 		{
    // 		  master_fprintf(file,"\n # %s%s\n\n",prop_name[contr_map[icontr].first],prop_name[contr_map[icontr].second]);
    // 		  for(int t=0;t<glbSize[0];t++)
    // 		    {
    // 		      int i=t+glbSize[0]*(icontr+contr_map.size()*(iflav+nflavs*jflav));
    // 		      master_fprintf(file, "%+16.16lg %+16.16lg\n",glb_contr[i][RE]*norm,glb_contr[i][IM]*norm);
    // 		    }
    // 		}
    // 	    }
    //   }
    
    // // vector_reset(ori_source[EVN]);
    // // vector_reset(ori_source[ODD]);
    // // for(int icol=0;icol<3;icol++) ori_source[EVN][0][icol][RE]=1;
    // // for(int par=0;par<2;par++) vector_copy(temp_source[par],ori_source[par]);
    // // MINV(M[0],0,temp_source);
    // // vector_reset(temp_source[EVN]);
    // // vector_reset(temp_source[ODD]);
    // // for(int par=0;par<2;par++)
    // //   {
    // // 	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
    // // 	  {
    // // 	    int ivol=loclx_of_loceo[par][ieo];
    // // 	    if(glb_coord_of_loclx[ivol][0]==glb_size[0]/2)
    // // 	      color_copy(temp_source[par][ieo],M[0][par][ieo]);
    // // 	  }
    // //  NISSA_PARALLEL_LOOP_END;
    // // 	set_borders_invalid(temp_source[par]);
    // //   }
    
    // // put_stag_phases(temp_source,form_stag_op_pattern(15,15));
    // // put_stag_phases(M[0],form_stag_op_pattern(15,15));
    // // MINV(M[1],0,temp_source);
    
    // // add_backfield_to_conf(conf,theory_pars.backfield[0]);
    // // vector_reset(loc_contr);
    // // for(int par=0;par<2;par++)
    // //   NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
    // // 	{
    // // 	  int ivol=loclx_of_loceo[par][ieo];
    // // 	  color f;
    // // 	  unsafe_su3_prod_color(f,conf[par][ieo][0],M[1][!par][loceo_neighup[par][ieo][0]]);
    // // 	  color b;
    // // 	  unsafe_su3_dag_prod_color(b,conf[par][ieo][0],M[1][!par][loceo_neighup[par][ieo][0]]);
	  
    // // 	  complex temp;
    // // 	  color_scalar_prod(temp,M[0][par][ieo],f);
    // // 	  complex_summassign(loc_contr[glb_coord_of_loclx[ivol][0]],temp);
    // // 	  color_scalar_prod(temp,M[0][par][ieo],b);
    // // 	  complex_summassign(loc_contr[glb_coord_of_loclx[ivol][0]],temp);
    // // 	}
    // // NISSA_PARALLEL_LOOP_END;
    // // rem_backfield_from_conf(conf,theory_pars.backfield[0]);
    
    // // for(int t=0;t<glb_size[0];t++)
    // //   master_printf("%d %lg %lg\n",t,loc_contr[t][RE],loc_contr[t][IM]);
    
    // //free
    // for(int par=0;par<2;par++)
    //   {
    // 	nissa_free(photon_field[par]);
    // 	nissa_free(temp_source[par]);
    // 	nissa_free(ori_source[par]);
    //   }
    // for(int i=0;i<nprop_t*nflavs;i++)
    //   for(int par=0;par<2;par++)
    // 	nissa_free(M[i][par]);
    // nissa_free(glb_contr);
    
    // close_file(file);
  }
}
