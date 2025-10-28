#include <math.h>

#include "nissa.h"

int main(int narg,char **arg)
{  initNissa();
  
  if(rank==0)
    {
      //check commutation
      for(int ig=1;ig<5;ig++)
	for(int jg=ig;jg<5;jg++)
	{
	  dirac_matr test1,test2;
	  dirac_prod(&test1,&(base_gamma[ig]),&(base_gamma[jg]));
	  dirac_prod(&test2,&(base_gamma[jg]),&(base_gamma[ig]));

	  dirac_matr test3;	
	  printf("%sx%s ",gtag[ig],gtag[jg]);
	  if(ig==jg)
	    {
	      dirac_subt(&test3,&test1,&test2);
	      printf("-");
	    }
	  else
	    {
	      dirac_summ(&test3,&test1,&test2);
	      printf("+");
	    }
	  printf(" %sx%s:\n",gtag[jg],gtag[ig]);

	  print_dirac(&test3);
	}
    }
    
  closeNissa();

  return 0;
}
