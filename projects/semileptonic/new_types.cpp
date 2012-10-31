#include <stdio.h>
#include <stdlib.h>
#include "../../src/nissa.h"

#include "./new_types.h"

char base_out_folder[1024];

// ######################################################### mass_res_group_t ###################################################

// ######################################################### theta_group_t ###################################################

// ######################################################### ape_smear_t ####################################################

// ########################################################### source_t #####################################################

void source_t::smear(gauge_conf_t &conf,gauss_smear_pars_t &pars)
{jacobi_smearing(eta,eta,conf.U,pars.kappa,pars.nterm,pars.coeff,pars.expnt);}

// ######################################################### prop_group_t ###################################################

void prop_group_t::get_ntheta_mass_r(int &ntheta,int &nmass,int &nr)
{
  ntheta=theta->ntheta;
  nmass=mass_res->nmass;
  nr=(which_r==R_BOTH)+1;
}

void prop_group_t::check_itheta_mass_r(int itheta,int imass,int r)
{
  int ntheta,nmass,nr;
  get_ntheta_mass_r(ntheta,nmass,nr);
  
  //check theta
  if(itheta<0||itheta>=ntheta) crash("asked for itheta=%d, prop group contains %d",itheta,ntheta);
  //check mass
  if(imass<0||imass>=nmass) crash("asked for imass=%d, prop group contains %d",imass,nmass);
  //check r
  if(r<0||r>=nr) crash("asked for r=%d, prop group contains %d",r,nr);
}

int prop_group_t::iprop(int itheta,int imass,int r)
{
  check_itheta_mass_r(itheta,imass,r);
  
  int ntheta,nmass,nr;
  get_ntheta_mass_r(ntheta,nmass,nr);
  
  return r+nr*(imass+nmass*itheta);
}

int prop_group_t::nprop()
{
  int ntheta,nmass,nr;
  get_ntheta_mass_r(ntheta,nmass,nr);

  return nr*ntheta*nmass;
}

void prop_group_t::create(theta_group_t &t,mass_res_group_t &m,TMR r)
{
  mass_res=&m;
  theta=&t;
  which_r=r;
  
  int n=nprop();
  S=nissa_malloc("S*",n,colorspinspin*);
  for(int i=0;i<n;i++) S[i]=nissa_malloc("S",loc_vol+bord_vol,colorspinspin);
}

void prop_group_t::read_pars(int ntheta_group,theta_group_t *theta_group,int nmass_res_group,mass_res_group_t *mass_res_group)
{
  //read the id of the list of masses and residues
  int imass_res_group;
  read_str_int("iMassResGroup",&imass_res_group);
  if(imass_res_group<0||imass_res_group>nmass_res_group)
    crash("mass_res group chosen %d is not in the interval [0,%d)",imass_res_group,nmass_res_group);
      
  //read the id of the list of theta
  int itheta_group;
  read_str_int("iThetaGroup",&itheta_group);
  if(itheta_group<0||itheta_group>ntheta_group)
    crash("theta group chosen %d is not in the interval [0,%d)",itheta_group,ntheta_group);
  
  //read which r
  int which_r;
  read_str_int("WhichR",&which_r);
  if(which_r<0||which_r>2)
    crash("r chosen %d not in the interval [0,2]",which_r);
  
  //allocate the propagators
  create(theta_group[itheta_group],mass_res_group[imass_res_group],(TMR)which_r);
}

void prop_group_t::get_inverting(source_t &source,gauge_conf_t &gauge_conf,int rotate_to_physical_basis)
{
  //get ntheta,nmass,nr
  int nmass,ntheta,nr;
  get_ntheta_mass_r(ntheta,nmass,nr);
  
  //allocate
  spincolor *temp_source=nissa_malloc("temp_source",loc_vol+bord_vol,spincolor);
  spincolor *temp_reco[2]={nissa_malloc("temp_reco",loc_vol+bord_vol,spincolor),nissa_malloc("temp_reco",loc_vol+bord_vol,spincolor)};
  spincolor *cgm_solution[nmass];
  for(int imass=0;imass<nmass;imass++)
    cgm_solution[imass]=nissa_malloc(combine("cgm_solution_%d",imass).c_str(),loc_vol+bord_vol,spincolor);
  
  for(int id=0;id<4;id++)
    { 
      //extract index of the source
      get_spincolor_from_colorspinspin(temp_source,source.eta,id);
      //put the g5
      unsafe_dirac_prod_spincolor(temp_source,base_gamma[5],temp_source);
      
      for(int itheta=0;itheta<ntheta;itheta++)
	{
	  //adapt the border condition
	  double th=theta->theta[itheta];
	  momentum_t mom={1,th,th,th};
	  gauge_conf.adapt_theta(mom);
                
	  //invert
	  int niter_max=100000;
	  inv_tmQ2_cgm(cgm_solution,gauge_conf.U,gauge_conf.kappa,mass_res->mass,nmass,niter_max,mass_res->residues,temp_source);
	  
	  for(int imass=0;imass<nmass;imass++)
	    {
	      //reconstruct the doublet
	      reconstruct_tm_doublet(temp_reco[0],temp_reco[1],gauge_conf.U,gauge_conf.kappa,mass_res->mass[imass],cgm_solution[imass]);
	      master_printf("Mass %d (%g) reconstructed \n",imass,mass_res->mass[imass]);
	      
	      //convert the id-th spincolor into the colorspinspin
	      for(int rdest=0;rdest<nr;rdest++)
		put_spincolor_into_colorspinspin(S[iprop(itheta,imass,rdest)],temp_reco[(nr==2)?rdest:which_r],id);
	    }
	}
    }
  
  //rotate if needed
  if(rotate_to_physical_basis)
    for(int itheta=0;itheta<ntheta;itheta++)
      for(int imass=0;imass<nmass;imass++)
	for(int rdest=0;rdest<nr;rdest++) //rotate opposite of D
	  rotate_vol_colorspinspin_to_physical_basis(S[iprop(itheta,imass,rdest)],!rdest,!rdest);
  
  //free
  for(int imass=0;imass<nmass;imass++) nissa_free(cgm_solution[imass]);
  nissa_free(temp_source);
  for(int r=0;r<2;r++) nissa_free(temp_reco[r]);
}

void prop_group_t::smear(gauge_conf_t &conf,gauss_smear_pars_t &pars)
{
  int n=nprop();
  
  for(int i=0;i<n;i++)
    jacobi_smearing(S[i],S[i],conf.U,pars.kappa,pars.nterm,pars.coeff,pars.expnt);
}

void prop_group_t::get_smearing(gauge_conf_t &conf,gauss_smear_pars_t &pars,prop_group_t &in)
{
  int nin=in.nprop();
  int nout=nprop();
  if(nout!=nin) crash("in group has %d elements, out prop has %d",nin,nout);
  
  for(int i=0;i<nin;i++)
    jacobi_smearing(S[i],in.S[i],conf.U,pars.kappa,pars.nterm,pars.coeff,pars.expnt);
}

void prop_group_t::get_reading(const char *ext_template_path,gauge_conf_t &conf,int load_reconstructing,int rotate_to_physical_basis)
{
  char template_path[1024];
  sprintf(template_path,"%s/%s",base_out_folder,ext_template_path);
  
  int ntheta,nmass,nr;
  get_ntheta_mass_r(ntheta,nmass,nr);

  for(int itheta=0;itheta<ntheta;itheta++)
    for(int imass=0;imass<nmass;imass++)
      {
	if(nr==2 && load_reconstructing)
	  {
	    double th=theta->theta[itheta];
	    momentum_t mom={1,th,th,th};
	    conf.adapt_theta(mom);
	    
	    int ip0=iprop(itheta,imass,0);
	    int ip1=iprop(itheta,imass,1);
	    colorspinspin *temp[2]={S[ip0],S[ip1]};
	    read_tm_colorspinspin_reconstructing(temp,combine(template_path,ip0).c_str(),NULL,conf.U,conf.kappa,mass_res->mass[imass]);
	  }
	else
	  for(int r=0;r<nr;r++)
	    {
	      int ip=iprop(itheta,imass,r);
	      read_colorspinspin(S[ip],combine(template_path,ip).c_str(),NULL);
	    }
	
	if(rotate_to_physical_basis)
	  for(int r=0;r<nr;r++)
	    {
	      int ip=iprop(itheta,imass,r);
	      rotate_vol_colorspinspin_to_physical_basis(S[ip],!r,!r);
	    }
      }
}

void prop_group_t::write(const char *ext_template_path,int save_reconstructing,int is_rotated,gauge_conf_t &gauge_conf)
{
  char template_path[1024];
  sprintf(template_path,"%s/%s",base_out_folder,ext_template_path);
  
  int ntheta,nmass,nr;
  get_ntheta_mass_r(ntheta,nmass,nr);

  for(int itheta=0;itheta<ntheta;itheta++)
    for(int imass=0;imass<nmass;imass++)
  
  for(int id=0;id<4;id++)
    {
      int ivol1=8,id1=2,ic1=1,ri1=1,mu1=1;
      int ip0=iprop(itheta,imass,0);
      int ip1=iprop(itheta,imass,1);
    }
  
  for(int itheta=0;itheta<ntheta;itheta++)
    for(int imass=0;imass<nmass;imass++)
      if(nr==2 && save_reconstructing)
	{
	  int ip0=iprop(itheta,imass,0);
	  int ip1=iprop(itheta,imass,1);
	  
	  double th=theta->theta[itheta];
	  momentum_t mom={1,th,th,th};
          gauge_conf.adapt_theta(mom);
	  master_printf("involved: %d %d\n",ip0,ip1);
	  write_tm_colorspinspin_anti_reconstructing(combine(template_path,ip0).c_str(),S[ip0],S[ip1],is_rotated,mass_res->mass[imass],64,gauge_conf.U,gauge_conf.kappa,gauge_conf.theta);
	}
      else
	for(int r=0;r<2;r++)
	  {
	    int ip=iprop(itheta,imass,r);
	    write_colorspinspin(combine(template_path,ip).c_str(),S[ip],64);
	  }  
}

// #################################### prop_group_command_t ##################################

void prop_group_command_t::read_command(prop_group_t &ext_prop_group_out,source_t *ext_source,prop_group_t *ext_prop_group_in,gauge_conf_t *ext_conf,gauss_smear_pars_t *ext_smear_pars)
{
  //copy reference group
  prop_group_out=&ext_prop_group_out;
  
  //read if to invert
  read_str_int("GetInverting",&get_inverting);
  
  if(get_inverting)
    {
      //we will not read
      get_reading=0;
      
      //read which source to use
      int obtain_inverting_isource;
      read_str_int("ObtainInvertingSourceId",&obtain_inverting_isource);
      source=ext_source+obtain_inverting_isource;
      
      //read_which conf to use
      int obtain_inverting_iconf;
      read_str_int("ObtainInvertingConfId",&obtain_inverting_iconf);
      conf=ext_conf+obtain_inverting_iconf;
      
      //read if to rotate
      read_str_int("RotateToPhysicalBasis",&rotate_to_physical_basis);
    }
  else
    {
      //read whether to load
      read_str_int("GetReading",&get_reading);
      
      if(get_reading)
	{
	  //read if to reconstruct
	  //put to default to avoid confusion
	  //read_str_int("ReadReconstructing",&load_reconstructing);
	  load_reconstructing=1;
	  
	  //read_which conf to use to reconstruct
	  if(load_reconstructing)
	    {
	      int reconstruct_using_iconf;
	      read_str_int("ReconstructUsingConfId",&reconstruct_using_iconf);
	      conf=ext_conf+reconstruct_using_iconf;
	    }
	  
	  //read if to rotate
	  read_str_int("RotateToPhysicalBasis",&rotate_to_physical_basis);
	}
      else
	{
	  //read in prop_group
	  int obtain_smearing_igroup;
	  read_str_int("ObtainSmearingGroupId",&obtain_smearing_igroup);
	  prop_group_in=ext_prop_group_in+obtain_smearing_igroup;
	  
	  //read operator to use for smearing
	  int obtain_smearing_with_ioperator;
	  read_str_int("ObtainSmearingWithOperatorId",&obtain_smearing_with_ioperator);
	  smear_pars=ext_smear_pars+obtain_smearing_with_ioperator;
	  
	  //read conf to be used to smear
	  int obtain_smearing_with_iconf;
	  read_str_int("ObtainSmearingWithConfId",&obtain_smearing_with_iconf);
	  conf=ext_conf+obtain_smearing_with_iconf;
	}
    }
  
  //if not to load, read whether to save
  if(!get_reading)
    {
      read_str_int("Write",&save);
      
      //putting to default to avoid confusion
      //if(save) read_str_int("WriteToBeReconstructed",&save_reconstructing);
      save_reconstructing=1;
    }

  //if to load or to write, read the template path
  if(get_reading||save) read_str_str("TemplatePath",template_path,1024);
}

//the conf and source are passed temporary
void prop_group_command_t::exec()
{
  if(get_reading) prop_group_out->get_reading(template_path,*conf,load_reconstructing,rotate_to_physical_basis);
  else
    {
      if(get_inverting) prop_group_out->get_inverting(*source,*conf,rotate_to_physical_basis);
      else prop_group_out->get_smearing(*conf,*smear_pars,*prop_group_in);
      //write
      if(save) prop_group_out->write(template_path,save_reconstructing,rotate_to_physical_basis,*conf);
    }
}

// ##################################### gauss_smear_pars_t ################################

void gauss_smear_pars_t::read()
{
  double ext_kappa;
  read_str_double("Kappa",&ext_kappa);
  set_kappa(ext_kappa);
  
  int ext_nterm;
  read_str_int("NTerms",&ext_nterm);
  create(ext_nterm);
  
  for(int iterm=0;iterm<ext_nterm;iterm++)
    {
      read_double(&coeff[iterm]);
      read_int(&expnt[iterm]);
      if(expnt[iterm]<0) crash("chosen exponent %d for %d term of smearing operator, select a positive one",expnt[iterm],iterm);
      if(iterm>0 && expnt[iterm]<expnt[iterm-1])
	crash("exponent %d (%d) smaller than %d (%d), please sort coefficient in ascending order",iterm,expnt[iterm],iterm-1,expnt[iterm-1]);
    }
}

// #################################### two_pts_corr_group_t ##############################

void two_pts_corr_group_t::add_corr(std::vector<two_pts_contr_pars_t> &buf_contr_list,std::vector<std::string> &buf_corr_name,const char *what)
{
  //use the same order as in the correlator (sink, source) ALWAYS
  
  if(strcmp("S0S0",what)==0)
    {
      buf_corr_name.push_back(std::string(what));
      buf_contr_list.push_back(two_pts_contr_pars_t(1,0,0,1.0));
    }
  if(strcmp("S0P5",what)==0)
    {
      buf_corr_name.push_back(std::string(what));
      buf_contr_list.push_back(two_pts_contr_pars_t(1,0,5,1.0));
    }
  if(strcmp("P5S0",what)==0)
    {
      buf_corr_name.push_back(std::string(what));
      buf_contr_list.push_back(two_pts_contr_pars_t(1,5,0,1.0));
    }
  if(strcmp("VKVK",what)==0)
    {
      buf_corr_name.push_back(std::string(what));
      buf_contr_list.push_back(two_pts_contr_pars_t(1,1,1,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,2,2,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,3,3,1.0/3));
    }
  if(strcmp("AKVK",what)==0)
    {
      buf_corr_name.push_back(std::string(what));
      buf_contr_list.push_back(two_pts_contr_pars_t(1,6,1,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,7,2,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,8,3,1.0/3));
    }
  if(strcmp("VKAK",what)==0)
    {
      buf_corr_name.push_back(std::string(what));
      buf_contr_list.push_back(two_pts_contr_pars_t(1,1,6,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,2,7,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,3,8,1.0/3));
    }
  if(strcmp("V0V0",what)==0)
    {
      buf_corr_name.push_back(std::string(what));
      buf_contr_list.push_back(two_pts_contr_pars_t(1,4,4,1.0));
    }
  if(strcmp("P5P5",what)==0)
    {
      buf_contr_list.push_back(two_pts_contr_pars_t(1,5,5,1));
      buf_corr_name.push_back(std::string(what));
    }
  if(strcmp("AKAK",what)==0)
    {
      buf_corr_name.push_back(std::string(what));
      buf_contr_list.push_back(two_pts_contr_pars_t(1,6,6,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,7,7,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,8,8,1.0/3));
    }
  if(strcmp("A0A0",what)==0)
    {
      buf_corr_name.push_back(std::string(what));
      buf_contr_list.push_back(two_pts_contr_pars_t(1,9,9,1.0));
    }
  if(strcmp("TKTK",what)==0)
    {
      buf_corr_name.push_back(std::string(what));
      buf_contr_list.push_back(two_pts_contr_pars_t(1,10,10,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,11,11,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,12,12,1.0/3));
    }
  if(strcmp("BKBK",what)==0)
    {
      buf_corr_name.push_back(std::string(what));
      buf_contr_list.push_back(two_pts_contr_pars_t(1,13,13,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,14,14,1.0/3));
      buf_contr_list.push_back(two_pts_contr_pars_t(0,15,15,1.0/3));
    }
}

void two_pts_corr_group_t::create(std::vector<two_pts_contr_pars_t> &buf_contr_list,std::vector<std::string> &buf_corr_name)
{
  //unbuffer corr
  corr_name=nissa_malloc("CorrList*",buf_corr_name.size(),char*);
  for(int icorr=0;icorr<buf_corr_name.size();icorr++)
    {
      corr_name[icorr]=(char*)malloc(10);
      sprintf(corr_name[icorr],"%s",buf_corr_name[icorr].c_str());
    }
  
  //unbuffer contr
  contr_list=nissa_malloc("ContrList",buf_contr_list.size(),two_pts_contr_pars_t);
  ncontr=buf_contr_list.size();
  for(int icontr=0;icontr<buf_contr_list.size();icontr++)
    contr_list[icontr]=buf_contr_list[icontr];  
}

void two_pts_corr_group_t::read()
{
  std::vector<two_pts_contr_pars_t> buf_contr_list;
  std::vector<std::string> buf_corr_name;
  
  //read the number of corr
  read_str_int("NCorr",&ncorr);
  
  //read the corr one by one
  for(int icorr=0;icorr<ncorr;icorr++)
    {
      char corr_name[1024];
      read_str(corr_name,1024);
      add_corr(buf_contr_list,buf_corr_name,corr_name);
    }
  
  create(buf_contr_list,buf_corr_name);
}

// ###################################### corr_command_t ###################################

void corr_command_t::read_prop_group_pair(int nprop_group,prop_group_t *prop_group,int ipair)
{
  //the first is reverted
  int first,second;
  read_int(&first);
  read_int(&second);
  
  if(first<0||first>=nprop_group) crash("first group %d must be in the interval [0,%d)",first,nprop_group);
  if(second<0||second>=nprop_group) crash("second group %d must be in the interval [0,%d)",second,nprop_group);
  
  pair_list[ipair]=prop_group_pair_t(prop_group[first],prop_group[second]);
}

void corr_command_t::read_corr_group(int ntwo_pts_corr_group_avail,two_pts_corr_group_t *ext_two_pts_corr_group)
{
  int two_pts_corr_igroup;
  read_str_int("TwoPtsCorrGroupId",&two_pts_corr_igroup);
  
  if(two_pts_corr_igroup<0||two_pts_corr_igroup>=ntwo_pts_corr_group_avail) crash("two pts corr group id %d must be in the range [0,%d)",two_pts_corr_igroup,ntwo_pts_corr_group_avail);
  two_pts_corr_group=ext_two_pts_corr_group+two_pts_corr_igroup;
}

void corr_command_t::read(int ntwo_pts_group_avail,two_pts_corr_group_t *ext_two_pts_corr_group,int nprop_group,prop_group_t *prop_group)
{
  read_str_str("OutputPath",path,1024);
  
  read_corr_group(ntwo_pts_group_avail,ext_two_pts_corr_group);
  
  read_str_int("NPropGroupPair",&nprop_group_pair);
  pair_list=nissa_malloc("PairList",nprop_group_pair,prop_group_pair_t);
  
  for(int ipair=0;ipair<nprop_group_pair;ipair++)
    read_prop_group_pair(nprop_group,prop_group,ipair);
}

void corr_command_t::exec()
{
  FILE *fout=open_file(combine("%s/%s",base_out_folder,path).c_str(),"w");
  
  for(int ipair=0;ipair<nprop_group_pair;ipair++)
    {
      master_printf("Starting contraction of group %d/%d\n",ipair,nprop_group_pair);
      
      int ntheta1=pair_list[ipair].first->theta->ntheta;
      double *theta1=pair_list[ipair].first->theta->theta;
      int ntheta2=pair_list[ipair].second->theta->ntheta;
      double *theta2=pair_list[ipair].second->theta->theta;
      int nmass1=pair_list[ipair].first->mass_res->nmass;
      double *mass1=pair_list[ipair].first->mass_res->mass;
      double *res1=pair_list[ipair].first->mass_res->residues;
      int nmass2=pair_list[ipair].second->mass_res->nmass;
      double *mass2=pair_list[ipair].second->mass_res->mass;
      double *res2=pair_list[ipair].second->mass_res->residues;
      
      int ncontr=two_pts_corr_group->ncontr;
      int ncorr=two_pts_corr_group->ncorr;
		
      //prepare the list of contractions
      int source_op[ncontr];
      int sink_op[ncontr];
      double coeff[ncontr];
      for(int icontr=0;icontr<ncontr;icontr++)
	{
	  source_op[icontr]=two_pts_corr_group->contr_list[icontr].source_op;
	  sink_op[icontr]=two_pts_corr_group->contr_list[icontr].sink_op;
	  coeff[icontr]=two_pts_corr_group->contr_list[icontr].coeff;
	}
      complex *buf=nissa_malloc("buf",2*ncontr*glb_size[0],complex);
      
      for(int itheta1=0;itheta1<ntheta1;itheta1++)
	for(int imass1=0;imass1<nmass1;imass1++)
	  for(int itheta2=0;itheta2<ntheta2;itheta2++)
	    for(int imass2=0;imass2<nmass2;imass2++)
	      {
		master_fprintf(fout," # group_pair=%d, m1=%lg th1=%lg res1=%lg (reverted), m2=%lg th2=%lg res2=%lg\n\n",
			       ipair,mass1[imass1],theta1[itheta1],res1[imass1],mass2[imass2],theta2[itheta2],res2[imass2]);
		
		//contract
		for(int r=0;r<2;r++)
		  {
		    int iprop1=pair_list[ipair].first->iprop(itheta1,imass1,r);
		    int iprop2=pair_list[ipair].second->iprop(itheta2,imass2,r);
		    meson_two_points_Wilson_prop(buf+r*glb_size[0]*ncontr,sink_op,pair_list[ipair].first->S[iprop1],source_op,pair_list[ipair].second->S[iprop2],ncontr);
		  }
		
		//add the contraction to build correlation functions
		int icontr=0,icorr=0;
		two_pts_contr_pars_t *contr=two_pts_corr_group->contr_list;
		char **corr_name=two_pts_corr_group->corr_name;
		do
		  {
		    //reset the corr
		    complex data[glb_size[0]];
		    memset(data,0,sizeof(complex)*glb_size[0]);
		    
		    //loop on contr
		    do
		      {
			for(int r=0;r<2;r++)
			  for(int t=0;t<glb_size[0];t++)
			    complex_summ_the_prod_double(data[t],buf[t+glb_size[0]*(icontr+r*ncontr)],0.5*coeff[icontr]);
			icontr++;
		      }
		    while(icontr!=ncontr && !contr->starting);
		    
		    master_fprintf(fout," # %s\n",corr_name[icorr]);
		    print_contraction_to_file(fout,-1,-1,data,shift,"",1);
		    master_fprintf(fout,"\n");
		    
		    icorr++;
		  }
		while(icorr!=ncorr);
	      }
     
      nissa_free(buf);
    }
  
  if(rank==0) fclose(fout);
}

// ####################################### gauge_conf_t ###################################

void gauge_conf_t::copy(gauge_conf_t &in)
{
  destroy();
  create();
  
  for(int mu=0;mu<4;mu++) theta[mu]=in.theta[mu];
  beta=in.beta;
  kappa=in.kappa;
  
  if(in.is_allocated()) vector_copy(U,in.U);
  else crash("copying from an unallocated conf");
}

void gauge_conf_t::read(const char *path)
{
  if(!is_allocated()) create();
  read_ildg_gauge_conf(U,path);
  reset_theta();
  
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(U));  
}

void gauge_conf_t::ape_smear(ape_smear_pars_t &ape_smear_pars)
{
  ape_spatial_smear_conf(U,U,ape_smear_pars.alpha,ape_smear_pars.niter);
  master_printf("smerded plaq: %.18g\n",global_plaquette_lx_conf(U));
}
