#include "nissa.hpp"

using namespace nissa;

static const int remez_extra_den_fact[2]={2,-1};

struct rational_approx_db_t
{
  static const int version=1;
  int npseudof_max;
  int nflav_max;
  int npoles_max;
  rat_approx_t *base;
  int ntot_approx(){return 2*npseudof_max*nflav_max*npoles_max;}
  rational_approx_db_t(int npseudof_max,int nflav_max,int npoles_max) : npseudof_max(npseudof_max),nflav_max(nflav_max),npoles_max(npoles_max)
  {    
    base=nissa_malloc("base",ntot_approx(),rat_approx_t);
  }
  ~rational_approx_db_t(){nissa_free(base);}
  rat_approx_t* get_appr(int npf,int nflav,int npol,int pfact){return base+pfact+2*(npol+npoles_max*(nflav+nflav_max*npf));}
};

int main(int narg,char **arg)
{
  int version=1;
  
  init_nissa(narg,arg);
  if(narg<2) crash("Use: %s input",arg[0]);
  
  //read input
  open_input(arg[1]);
  int npseudof_max,nflav_max,npoles_max;
  read_str_int("NPseudofMax",&npseudof_max);
  read_str_int("NFlavMax",&nflav_max);
  read_str_int("NPolesMax",&npoles_max);
  char outpath[1024];
  read_str_str("Outpath",outpath,1024);
  close_input();
  
  rational_approx_db_t rat_appr_db(npseudof_max,nflav_max,npoles_max);

  for(int npf=0;npf<npseudof_max;npf++)
    for(int nfl=0;nfl<nflav_max;nfl++)
      for(int npol=0;npol<npoles_max;npol++)
	for(int pfact=0;pfact<2;pfact++)
	  {
	    //allocate
	    rat_approx_t *rat=rat_appr_db.get_appr(npf,nfl,npol,pfact);
	    rat_approx_create(rat,npol,"");
	    snprintf(rat->name,20,"x^(%d/%+d)[%d]",nfl,remez_extra_den_fact[pfact]*npf,npol);
	    
	    //set pars
	    rat->num=nfl;
	    rat->den=remez_extra_den_fact[pfact]*npf;
	    
	  }
  
  
  //open the output and write this infos
  FILE *fout=open_file(outpath,"w");
  
  fclose(fout);
  
  close_nissa();
  
  return 0;
}
