#include <string.h>

#include "new_types_definitions.h"
#include "complex.h"
#include "spin.h"

//saturate two anti-simmetric tensors
void as2t_saturate(complex out,as2t a,as2t b)
{
  unsafe_complex_prod(out,a[0],b[0]);
  for(int munu=1;munu<6;munu++) complex_summ_the_prod(out,a[munu],b[munu]);
}

//Print a spinspin
void print_spinspin(spinspin s)
{
  for(int id1=0;id1<4;id1++)
    {
      for(int id2=0;id2<4;id2++) printf("%+016.16f,%+016.16f\t",s[id1][id2][0],s[id1][id2][1]);
      printf("\n");
    }
}

//trace of the product with a dirac matr of a spinspin
void summ_the_trace_prod_dirac_spinspin(complex c,dirac_matr *a,spinspin b)
{
  for(int id=0;id<4;id++)
    complex_summ_the_prod(c,a->entr[id],b[a->pos[id]][id]);
}
void trace_prod_dirac_spinspin(complex c,dirac_matr *a,spinspin b)
{
  c[0]=c[1]=0;
  summ_the_trace_prod_dirac_spinspin(c,a,b);
}

//prouduct of two spinspins
void spinspin_spinspindag_prod(spinspin out,spinspin a,spinspin b)
{
  //This is the line on the matrix
  memset(out,0,sizeof(spinspin));
  for(int id1=0;id1<4;id1++)
    for(int id2=0;id2<4;id2++)
      for(int id=0;id<4;id++)
	complex_summ_the_conj2_prod(out[id1][id2],a[id1][id],b[id2][id]);
}


//Get a spincolor from a colorspinspin
//In a spinspin the sink index runs slower than the source
void get_spincolor_from_colorspinspin(spincolor out,colorspinspin in,int id_source)
{
  for(int ic_sink=0;ic_sink<3;ic_sink++)
    for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
      {
	out[id_sink][ic_sink][0]=in[ic_sink][id_sink][id_source][0];
	out[id_sink][ic_sink][1]=in[ic_sink][id_sink][id_source][1];
      }
}

//Get a spincolor from a su3spinspin
void get_spincolor_from_su3spinspin(spincolor out,su3spinspin in,int id_source,int ic_source)
{
  for(int ic_sink=0;ic_sink<3;ic_sink++)
    for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
      {
	out[id_sink][ic_sink][0]=in[ic_sink][ic_source][id_sink][id_source][0];
	out[id_sink][ic_sink][1]=in[ic_sink][ic_source][id_sink][id_source][1];
      }
}

//Put a spincolor into a colorspinspin
void put_spincolor_into_colorspinspin(colorspinspin out,spincolor in,int id_source)
{
  for(int ic_sink=0;ic_sink<3;ic_sink++)
    for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
      {
	out[ic_sink][id_sink][id_source][0]=in[id_sink][ic_sink][0];
	out[ic_sink][id_sink][id_source][1]=in[id_sink][ic_sink][1];
      }
}

//Put a spincolor into a su3spinspin
void put_spincolor_into_su3spinspin(su3spinspin out,spincolor in,int id_source,int ic_source)
{
  for(int ic_sink=0;ic_sink<3;ic_sink++)
    for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
      {
	out[ic_sink][ic_source][id_sink][id_source][0]=in[id_sink][ic_sink][0];
	out[ic_sink][ic_source][id_sink][id_source][1]=in[id_sink][ic_sink][1];
      }
}
