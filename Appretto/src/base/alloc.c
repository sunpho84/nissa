#pragma once

//allocate vectors of the required length
char *allocate_vector(int length,char *tag)
{
  char *out=(char*)malloc(length);
  if(out==NULL && rank==0)
    {
      fprintf(stderr,"Error during allocation of %s\n",tag);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  return out;
}

spincolor *allocate_spincolor(int length,char *tag){return (spincolor*)allocate_vector(length*sizeof(spincolor),tag);}
redspincolor *allocate_redspincolor(int length,char *tag){return (redspincolor*)allocate_vector(length*sizeof(spincolor),tag);}
quad_su3 *allocate_quad_su3(int length,char *tag){return (quad_su3*)allocate_vector(length*sizeof(quad_su3),tag);}
su3 *allocate_su3(int length,char *tag){return (su3*)allocate_vector(length*sizeof(su3),tag);}
as2t_su3 *allocate_as2t_su3(int length,char *tag){return (as2t_su3*)allocate_vector(length*sizeof(as2t_su3),tag);}
colorspinspin *allocate_colorspinspin(int length,char *tag){return (colorspinspin*)allocate_vector(length*sizeof(colorspinspin),tag);}
su3spinspin *allocate_su3spinspin(int length,char *tag){return (su3spinspin*)allocate_vector(length*sizeof(su3spinspin),tag);}
