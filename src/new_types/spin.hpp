#ifndef _SPIN_HPP
#define _SPIN_HPP

#include <string.h>
#include <stdio.h>

#include "complex.hpp"
#include "new_types_definitions.hpp"
#include "spin.hpp"

#include "routines/math_routines.hpp"

namespace nissa
{
  //Print a spin
  inline void spin_print(spin s)
  {
    for(int id=0;id<4;id++) printf("%+016.16le,%+016.16le\t",s[id][0],s[id][1]);
    printf("\n");
  }
  
  inline void spin_put_to_zero(spin s)
  {for(int i=0;i<4;i++) complex_put_to_zero(s[i]);}
  
  inline void spin_copy(spin out,spin in)
  {for(int i=0;i<4;i++) complex_copy(out[i],in[i]);}
  inline void spin_conj(spin out,spin in)
  {for(int i=0;i<4;i++) complex_conj(out[i],in[i]);}
  
  inline void spin_summ(spin a,spin b,spin c)
  {for(int i=0;i<4;i++) complex_summ(a[i],b[i],c[i]);}
  inline void spin_summassign(spin a,spin b)
  {spin_summ(a,a,b);}
  inline void spin_subt(spin a,spin b,spin c)
  {for(int i=0;i<4;i++) complex_subt(a[i],b[i],c[i]);}
  inline void spin_subtassign(spin a,spin b)
  {spin_subt(a,a,b);}
  
  inline void spin_prod_double(spin a,spin b,double c)
  {for(int i=0;i<4;i++) complex_prod_double(a[i],b[i],c);}
  inline void spin_prodassign_double(spin a,double b)
  {spin_prod_double(a,a,b);}
  
  inline void spin_summ_the_complex_prod(spin a,spin b,complex c)
  {for(int i=0;i<4;i++) complex_summ_the_prod(a[i],b[i],c);}
  inline void spin_subt_the_complex_prod(spin a,spin b,complex c)
  {for(int i=0;i<4;i++) complex_subt_the_prod(a[i],b[i],c);}
  
  inline void spin_summ_the_complex_conj2_prod(spin a,spin b,complex c)
  {for(int i=0;i<4;i++) complex_summ_the_conj2_prod(a[i],b[i],c);}
  inline void spin_subt_the_complex_conj2_prod(spin a,spin b,complex c)
  {for(int i=0;i<4;i++) complex_subt_the_conj2_prod(a[i],b[i],c);}
  
  inline void spinspin_copy(spinspin b,spinspin a) {for(int i=0;i<4;i++) spin_copy(b[i],a[i]);}
  
  inline void unsafe_spinspin_hermitian(spinspin b,spinspin a)
  {
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	complex_conj(b[id1][id2],a[id2][id1]);
  }
  inline void safe_spinspin_hermitian(spinspin b,spinspin a)
  {
    spinspin temp;
    unsafe_spinspin_hermitian(temp,a);
    spinspin_copy(b,temp);
  }
  
  inline void spinspin_put_to_zero(spinspin a)
  {for(int i=0;i<4;i++) spin_put_to_zero(a[i]);}
  
  inline void spinspin_put_to_diag(spinspin a,complex b)
  {
    spinspin_put_to_zero(a);
    for(int id=0;id<4;id++)
      complex_copy(a[id][id],b);
  }
  inline void spinspin_put_to_diag(spinspin a,double b)
  {complex c={b,0};spinspin_put_to_diag(a,c);}
  inline void spinspin_put_to_idiag(spinspin a,double b)
  {complex c={0,b};spinspin_put_to_diag(a,c);}

  inline void spin_direct_prod(spinspin out,spin a,spin b)
  {
    for(int id=0;id<4;id++)
      for(int jd=0;jd<4;jd++)
	unsafe_complex_prod(out[id][jd],a[id],b[jd]);
  }
  
  inline void spinspin_put_to_id(spinspin a)
  {
    complex o={1,0};
    spinspin_put_to_diag(a,o);
  }
  
  inline void spinspin_summ(spinspin a,spinspin b,spinspin c)
  {for(int id1=0;id1<4;id1++) for(int id2=0;id2<4;id2++) complex_summ(a[id1][id2],b[id1][id2],c[id1][id2]);}
  inline void spinspin_subt(spinspin a,spinspin b,spinspin c)
  {for(int id1=0;id1<4;id1++) for(int id2=0;id2<4;id2++) complex_subt(a[id1][id2],b[id1][id2],c[id1][id2]);}
  inline void spinspin_summassign(spinspin a,spinspin b)
  {spinspin_summ(a,a,b);}
  inline void spinspin_subtassign(spinspin a,spinspin b)
  {spinspin_subt(a,a,b);}
  
  inline void spinspin_prod_double(spinspin a,spinspin b,double c)
  {for(int id1=0;id1<4;id1++) for(int id2=0;id2<4;id2++) complex_prod_double(a[id1][id2],b[id1][id2],c);}
  inline void spinspin_prodassign_double(spinspin a,double b)
  {spinspin_prod_double(a,a,b);}
  inline void spinspin_summ_the_prod_double(spinspin a,spinspin b,double c)
  {for(int id1=0;id1<4;id1++) for(int id2=0;id2<4;id2++) complex_summ_the_prod_double(a[id1][id2],b[id1][id2],c);}
  inline void spinspin_prod_idouble(spinspin a,spinspin b,double c)
  {for(int id1=0;id1<4;id1++) for(int id2=0;id2<4;id2++) complex_prod_idouble(a[id1][id2],b[id1][id2],c);}
  inline void spinspin_prodassign_idouble(spinspin a,double b)
  {spinspin_prod_idouble(a,a,b);}
  inline void spinspin_summ_the_prod_idouble(spinspin a,spinspin b,double c)
  {for(int id1=0;id1<4;id1++) for(int id2=0;id2<4;id2++) complex_summ_the_prod_idouble(a[id1][id2],b[id1][id2],c);}
  
  inline void spinspin_summ_the_complex_prod(spinspin a,spinspin b,complex c)
  {for(int id1=0;id1<4;id1++) for(int id2=0;id2<4;id2++) complex_summ_the_prod(a[id1][id2],b[id1][id2],c);}
  inline void spinspin_subt_the_complex_prod(spinspin a,spinspin b,complex c)
  {for(int id1=0;id1<4;id1++) for(int id2=0;id2<4;id2++) complex_subt_the_prod(a[id1][id2],b[id1][id2],c);}
  inline void spinspin_summassign_the_complex_prod(spinspin a,complex b){spinspin_summ_the_complex_prod(a,a,b);}
  inline void spinspin_subtassign_the_complex_prod(spinspin a,complex b){spinspin_subt_the_complex_prod(a,a,b);}
  
  inline void spinspin_summ_the_complex_conj2_prod(spinspin a,spinspin b,complex c)
  {for(int id1=0;id1<4;id1++) for(int id2=0;id2<4;id2++) complex_summ_the_conj2_prod(a[id1][id2],b[id1][id2],c);}
  inline void spinspin_subt_the_complex_conj2_prod(spinspin a,spinspin b,complex c)
  {for(int id1=0;id1<4;id1++) for(int id2=0;id2<4;id2++) complex_subt_the_conj2_prod(a[id1][id2],b[id1][id2],c);}
  
  inline void unsafe_spinspin_prod_complex(spinspin a,spinspin b,complex c)
  {
    spinspin_put_to_zero(a);
    spinspin_summ_the_complex_prod(a,b,c);
  }
  inline void safe_spinspin_prod_complex(spinspin a,spinspin b,complex c)
  {
    spinspin d;
    unsafe_spinspin_prod_complex(d,b,c);
    spinspin_copy(a,d);
  }

  inline void unsafe_spinspin_prod_complex_conj2(spinspin a,spinspin b,complex c)
  {
    complex d;
    complex_conj(d,c);
    unsafe_spinspin_prod_complex(a,b,d);
  }
  inline void safe_spinspin_prod_complex_conj2(spinspin a,spinspin b,complex c)
  {
    complex d;
    complex_conj(d,c);
    safe_spinspin_prod_complex(a,b,d);
  }
  
  //saturate two anti-simmetric tensors
  inline void as2t_saturate(complex out,as2t a,as2t b)
  {
    unsafe_complex_prod(out,a[0],b[0]);
    for(int munu=1;munu<6;munu++) complex_summ_the_prod(out,a[munu],b[munu]);
  }
  
  //Print a spinspin
  inline void spinspin_print(spinspin s)
  {
    for(int id1=0;id1<4;id1++)
      {
	for(int id2=0;id2<4;id2++) printf("%+016.16le,%+016.16le\t",s[id1][id2][0],s[id1][id2][1]);
	printf("\n");
      }
  }
  
  //trace of the product with a dirac matr of a spinspin
  inline void summ_the_trace_dirac_prod_spinspin(complex c,dirac_matr *a,spinspin b)
  {
    for(int id=0;id<4;id++)
      complex_summ_the_prod(c,a->entr[id],b[a->pos[id]][id]);
  }
  inline void trace_dirac_prod_spinspin(complex c,dirac_matr *a,spinspin b)
  {
    c[0]=c[1]=0;
    summ_the_trace_dirac_prod_spinspin(c,a,b);
  }
  
  inline void summ_the_trace_spinspin(complex c,spinspin a)
  {for(int id=0;id<4;id++) complex_summassign(c,a[id][id]);}
  inline void trace_spinspin(complex c,spinspin a)
  {
    c[0]=c[1]=0;
    summ_the_trace_spinspin(c,a);
  }
  
  //product of two spinspins
  inline void spinspin_summ_the_spinspin_dag_prod(spinspin out,spinspin a,spinspin b)
  {
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	for(int id=0;id<4;id++)
	  complex_summ_the_conj2_prod(out[id1][id2],a[id1][id],b[id2][id]);
  }
  inline void unsafe_spinspin_prod_spinspin_dag(spinspin out,spinspin a,spinspin b)
  {
    memset(out,0,sizeof(spinspin));
    spinspin_summ_the_spinspin_dag_prod(out,a,b);
  }
  inline void safe_spinspin_prod_spinspin_dag(spinspin out,spinspin a,spinspin b)
  {
    spinspin c;
    unsafe_spinspin_prod_spinspin_dag(c,a,b);
    memcpy(out,c,sizeof(spinspin));
  }
  inline void spinspin_summ_the_spinspin_prod(spinspin out,spinspin a,spinspin b)
  {
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	for(int id=0;id<4;id++)
	  complex_summ_the_prod(out[id1][id2],a[id1][id],b[id][id2]);
  }
  inline void spinspin_subt_the_spinspin_prod(spinspin out,spinspin a,spinspin b)
  {
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	for(int id=0;id<4;id++)
	  complex_subt_the_prod(out[id1][id2],a[id1][id],b[id][id2]);
  }
  inline void unsafe_spinspin_prod_spinspin(spinspin out,spinspin a,spinspin b)
  {
    memset(out,0,sizeof(spinspin));
    spinspin_summ_the_spinspin_prod(out,a,b);
  }
  inline void safe_spinspin_prod_spinspin(spinspin out,spinspin a,spinspin b)
  {
    spinspin c;
    unsafe_spinspin_prod_spinspin(c,a,b);
    memcpy(out,c,sizeof(spinspin));
  }
  
  inline double real_part_of_trace_spinspin_prod_spinspin_dag(spinspin a,spinspin b)
  {
    double t=0;
    
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	t+=a[id1][id2][0]*b[id1][id2][0]+a[id1][id2][1]*b[id1][id2][1];
    
    return t;
  }
  
  inline double spinspin_norm2(spinspin a)
  {return real_part_of_trace_spinspin_prod_spinspin_dag(a,a);}
  
  inline void trace_spinspin_prod_spinspin(complex c,spinspin a,spinspin b)
  {
    complex_put_to_zero(c);
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	complex_summ_the_prod(c,a[id1][id2],b[id2][id1]);
  }
  
  //prouduct of spinspin and spin
  inline void unsafe_spinspin_prod_spin(spin out,spinspin a,spin b)
  {
    memset(out,0,sizeof(spin));
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	complex_summ_the_prod(out[id1],a[id1][id2],b[id2]);
  }
  inline void safe_spinspin_prod_spin(spin out,spinspin a,spin b)
  {
    spin c;
    unsafe_spinspin_prod_spin(c,a,b);
    memcpy(out,c,sizeof(spin));
  }
  
  //prouduct of spin and spinspin
  inline void unsafe_spin_prod_spinspin(spin out,spin a,spinspin b)
  {
    memset(out,0,sizeof(spin));
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	complex_summ_the_prod(out[id1],a[id2],b[id2][id1]);
  }
  inline void safe_spin_prod_spinspin(spin out,spin a,spinspin b)
  {
    spin c;
    unsafe_spin_prod_spinspin(c,a,b);
    memcpy(out,c,sizeof(spin));
  }
  
  //Get a color from a colorspinspin
  inline void get_color_from_colorspinspin(color out,colorspinspin in,int id1,int id2)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[ic_sink],in[ic_sink][id1][id2]);}
  
  //Get a color from a spincolor
  inline void get_color_from_spincolor(color out,spincolor in,int id)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[ic_sink],in[id][ic_sink]);}
  
  //Get a spincolor from a colorspinspin
  //In a spinspin the sink index runs slower than the source
  inline void get_spincolor_from_colorspinspin(spincolor out,colorspinspin in,int id_source)
  {
    for(int ic_sink=0;ic_sink<NCOL;ic_sink++)
      for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
	complex_copy(out[id_sink][ic_sink],in[ic_sink][id_sink][id_source]);
  }
  
  //Get a spincolor from a su3spinspin
  inline void get_spincolor_from_su3spinspin(spincolor out,su3spinspin in,int id_source,int ic_source)
  {
    for(int ic_sink=0;ic_sink<NCOL;ic_sink++)
      for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
	complex_copy(out[id_sink][ic_sink],in[ic_sink][ic_source][id_sink][id_source]);
  }
  
  //Get a color from a su3
  inline void get_color_from_su3(color out,su3 in,int ic_source)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[ic_sink],in[ic_sink][ic_source]);}
  
  //Put a color into a colorspinspin
  inline void put_color_into_colorspinspin(colorspinspin out,color in,int id1,int id2)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[ic_sink][id1][id2],in[ic_sink]);}
  
  //Put a color into a spincolor
  inline void put_color_into_spincolor(spincolor out,color in,int id)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[id][ic_sink],in[ic_sink]);}
  
  //Put a spincolor into a colorspinspin
  inline void put_spincolor_into_colorspinspin(colorspinspin out,spincolor in,int id_source)
  {
    for(int ic_sink=0;ic_sink<NCOL;ic_sink++)
      for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
	complex_copy(out[ic_sink][id_sink][id_source],in[id_sink][ic_sink]);
  }
  
  //Put a spincolor into a su3spinspin
  inline void put_spincolor_into_su3spinspin(su3spinspin out,spincolor in,int id_source,int ic_source)
  {
    for(int ic_sink=0;ic_sink<NCOL;ic_sink++)
      for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
	complex_copy(out[ic_sink][ic_source][id_sink][id_source],in[id_sink][ic_sink]);
  }
  
  //Put a spincolor into a su3
  inline void put_color_into_su3(su3 out,color in,int ic_source)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[ic_sink][ic_source],in[ic_sink]);}
  
  //dirac*spinr
  inline void unsafe_dirac_prod_spin(spin out,dirac_matr *m,spin in)
  {for(int id1=0;id1<4;id1++) unsafe_complex_prod(out[id1],m->entr[id1],in[m->pos[id1]]);}
  inline void safe_dirac_prod_spin(spin out,dirac_matr *m,spin in){spin tmp;unsafe_dirac_prod_spin(tmp,m,in);spin_copy(out,tmp);}
  inline void unsafe_spin_prod_dirac(spin out,spin in,dirac_matr *m)
  {spin_put_to_zero(out);for(int id1=0;id1<4;id1++) complex_summ_the_prod(out[m->pos[id1]],in[id1],m->entr[id1]);}
  inline void safe_spin_prod_spin(spin out,spin in,dirac_matr *m){spin tmp;unsafe_spin_prod_dirac(tmp,in,m);spin_copy(out,tmp);}
  
  //Rotate left and right by (1+-ig5)/sqrt(2)
  //We distinguish for types of rotations, ir=rsink*2+rsource
  //We divide the spinspin into 4 blocks, according to id<2
  // -------------------------------
  // | div |  0  |  1  |  2  |  3  |
  // -------------------------------
  // | 0 1 | + 1 | 1 + | 1 - | - 1 |
  // | 2 3 | 1 - | - 1 | + 1 | 1 + |
  // -------------------------------
  // so just have to specify which is the block which rotate as +-i
  inline void rotate_spinspin_to_physical_basis(spinspin s,int rsi,int rso)
  {
    const int list_prb[4]={0,1,2,3},list_mrb[4]={3,2,1,0}; //plus and minus rotating blocks
    const int so_shft[4]={0,2,0,2},si_shft[4]={0,0,2,2};   //start of dirac indexes defining blocks
    
    int ir=rsi*2+rso,prb=list_prb[ir],mrb=list_mrb[ir];
    
    for(int dso=0;dso<2;dso++)
      for(int dsi=0;dsi<2;dsi++)
	{
	  int pso=dso+so_shft[prb],psi=dsi+si_shft[prb];
	  int mso=dso+so_shft[mrb],msi=dsi+si_shft[mrb];
	  
	  // rotation with +,-
	  assign_complex_prod_i(s[pso][psi]);
	  assign_complex_prod_minus_i(s[mso][msi]);
	}
  }
  
  //compute the determinant of spinspin matrix
  inline void spinspin_det(complex d,spinspin s)
  {matrix_determinant(d,(complex*)s,4);}
}

#endif
