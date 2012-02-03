#include <include.h>
#include <iostream>

using namespace std;

int nspec=2;
int ntheta=3;
int nmass=2;
int nmass_3pts[2]={2,1};
int T=48,L=24;
int njack=13;

int CSPEC=0;
int LSPEC=1;

int CHARM=1;
int LIGHT=0;

int STAND=0;
int POSTH=1;
int NEGTH=2;

int REAL=0;
int IMAG=1;

int R0=0;
int R1=1;

class info_3pts
{
public:
  int ispec;
  int ith_s1;
  int im_s1;
  int ith_s0;
  int im_s0;
  int reim;
};

class info_2pts
{
public:
  int ispec;
  int r_spec;
  int im_spec;
  int ith_seq;
  int r_seq;
  int im_seq;
  int reim;
 };

int icombo_3pts(int ispec,int ith_s0,int im_s0,int ith_s1,int im_s1,int reim)
{
  im_s1-=nmass-nmass_3pts[ispec];
  return reim+2*(im_s0+nmass*(ith_s0+ntheta*(im_s1+nmass_3pts[ispec]*ith_s1)));
}

int icombo_3pts(info_3pts &desc)
{return icombo_3pts(desc.ispec,desc.ith_s0,desc.im_s0,desc.ith_s1,desc.im_s1,desc.reim);}


jvec load_corr_3pts(const char *name,info_3pts &desc)
{
  int combo=icombo_3pts(desc);
  cout<<"Loading 3pts corr "<<combo<<" for spec "<<desc.ispec<<endl;
  return jvec_load(combine("3pts_sp%d_%s_30_30",desc.ispec,name).c_str(),T,njack,combo);
}


int icombo_2pts(int r_spec,int im_spec,int ith_seq,int r_seq,int im_seq,int reim)
{return reim+2*(r_spec+2*(im_spec+nmass*(r_seq+2*(im_seq+nmass*ith_seq))));}

int icombo_2pts(info_2pts &desc)
{return icombo_2pts(desc.r_spec,desc.im_spec,desc.ith_seq,desc.r_seq,desc.im_seq,desc.reim);}

jvec load_corr_2pt(const char *name,info_2pts &desc)
{
  int combo=icombo_2pts(desc);
  cout<<"Loading 2pts corr "<<combo<<endl;
  return jvec_load(combine("2pts_%s_30_30",name).c_str(),T,njack,combo);
}

void compute_ZV(info_3pts &three,info_2pts &two)
{
  jvec V0P5=load_corr_3pts("V0P5",three).simmetrized(-1);
  jvec P5P5=load_corr_2pt("P5P5",two);
  
  jack cenP5P5=P5P5[L]/2;
  
  {
    ofstream out("V0P5.xmg");
    out<<"@type xydy"<<endl;
    out<<V0P5<<endl;
  }
  
  {
    ofstream out("P5P5.xmg");
    out<<"@type xydy"<<endl;
    out<<effective_mass(P5P5.simmetrized(1))<<endl;
  }
  
  jack ZV=-cenP5P5/constant_fit(V0P5,8,16);
  cout<<ZV<<endl;
}

int main()
{
  ////////////////////////////////////////////////////////////////
  info_3pts STAND_D_STAND_D_CSPEC;
  /*
      / \
   l /   \ l
    /  c  \
    ------- */
  STAND_D_STAND_D_CSPEC.ispec=CSPEC;
  STAND_D_STAND_D_CSPEC.ith_s0=STAND;
  STAND_D_STAND_D_CSPEC.ith_s1=STAND;
  STAND_D_STAND_D_CSPEC.im_s0=LIGHT;
  STAND_D_STAND_D_CSPEC.im_s1=LIGHT;
  STAND_D_STAND_D_CSPEC.reim=REAL;
  
  ////////////////////////////////////////////////////////////////
  info_3pts POSTH_D_POSTH_D_CSPEC=STAND_D_STAND_D_CSPEC;
  /*
      / \
   l+/   \l+
    /  c  \
    ------- */
  POSTH_D_POSTH_D_CSPEC.ith_s0=POSTH;
  POSTH_D_POSTH_D_CSPEC.ith_s1=POSTH;
  
  ///////////////////////////////////////////////////////////////
  info_3pts STAND_ETA_STAND_ETA=STAND_D_STAND_D_CSPEC;
  /*
      / \
   c /   \ c
    /  c  \
    ------- */
  STAND_ETA_STAND_ETA.im_s0=CHARM;
  STAND_ETA_STAND_ETA.im_s1=CHARM;

  ///////////////////////////////////////////////////////////////
  info_3pts STAND_D_STAND_D_LSPEC=STAND_ETA_STAND_ETA;
  /*
      / \
   c /   \ c
    /  l  \
    ------- */
  STAND_D_STAND_D_LSPEC.ispec=LSPEC;
  
  ////////////////////////////////////////////////////////////////
  info_2pts STAND_D_CSPEC_00;
  /*
    l0 -------
    
    c0 ------- */
  STAND_D_CSPEC_00.ispec=CSPEC;
  STAND_D_CSPEC_00.r_spec=R0;
  STAND_D_CSPEC_00.im_spec=CHARM;
  STAND_D_CSPEC_00.ith_seq=STAND;
  STAND_D_CSPEC_00.r_seq=R0;
  STAND_D_CSPEC_00.im_seq=LIGHT;
  STAND_D_CSPEC_00.reim=REAL;
  
  ////////////////////////////////////////////////////////////////
  info_2pts STAND_D_LSPEC_00=STAND_D_CSPEC_00;
  /*
    c0 -------
    
    l0 ------- */
  STAND_D_LSPEC_00.im_spec=LIGHT;
  STAND_D_LSPEC_00.im_seq=CHARM;
  STAND_D_LSPEC_00.ispec=LSPEC;

  ////////////////////////////////////////////////////////////////
  info_2pts POSTH_D_CSPEC_00=STAND_D_CSPEC_00;
  /*
    l0 -------
    
    c0 ------- */
  POSTH_D_CSPEC_00.ith_seq=POSTH;
  
  ////////////////////////////////////////////////////////////////
  info_2pts STAND_ETA_00=STAND_D_CSPEC_00;
  /*
    c0 -------
    
    c0 ------- */
  STAND_ETA_00.im_seq=CHARM;
  
  ////////////////////////////////////////////////////////////////
  
  compute_ZV(POSTH_D_POSTH_D_CSPEC,POSTH_D_CSPEC_00);
  compute_ZV(STAND_D_STAND_D_CSPEC,STAND_D_CSPEC_00);
  compute_ZV(STAND_ETA_STAND_ETA,STAND_ETA_00);
  compute_ZV(STAND_D_STAND_D_LSPEC,STAND_D_LSPEC_00);
  
  return 0;
}
