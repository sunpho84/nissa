#include <nissa.hpp>

// #include "base/quda_bridge.hpp"

using namespace nissa;

// double rel_diff_norm(spincolor *test,spincolor *ref)
// {
//   double_vector_subtassign((double*)test,(double*)ref,locVol*sizeof(spincolor)/sizeof(double));
//   const double norm2_diff=double_vector_glb_norm2(test,locVol);
//   const double norm2_ref=double_vector_glb_norm2(ref,locVol);
//   const double res=sqrt(norm2_diff/norm2_ref);
  
//   return res;
// }

DECLARE_TRANSPOSABLE_COMP(Spin,int,NDIRAC,spin);
DECLARE_TRANSPOSABLE_COMP(Color,int,NCOL,color);

void in_main(int narg,char **arg)
{
  const int T=16,L=16;
  
  init_grid(T,L);
  
  LxField<quad_su3> conf("conf");
  
  {
    DynamicTens<OfComps<SpinRow,LocEoSite>,double, MemoryType::CPU> d(std::make_tuple(locEoSite(locVolh)));
    auto e=d.getReadable();
    // e(spinRow(1))=1.0;
  }
  
  {
    constexpr StackTens<OfComps<ComplId>,double> I{0,1};
    
    verbosity_lv=3;
    Field2<OfComps<SpinRow,ComplId>,double> df(WITH_HALO);
    df.updateHalo();
    auto fl=df.flatten();
    fl(decltype(df)::FlattenedInnerComp(0));
    Field2<OfComps<ComplId>,double> ds;
    // EoField2<OfComps<SpinRow>,double> df2;
    //auto rdf=df.getWritable();
    // auto rdf2=df2.getWritable();
    df=I;
    master_printf("written 0? %lg\n",df(locLxSite(0),spinRow(0),reIm(0)));
    master_printf("written 1? %lg\n",df(locLxSite(0),spinRow(0),reIm(1)));
    df=df+I;
    ds=trace(dag(df(spinRow(0))));
    master_printf("written 2? %lg\n",df(locLxSite(0),spinRow(0),reIm(1)));
    master_printf("copied -2? %lg\n",ds(locLxSite(0),reIm(1)));
    // rdf(spinRow(0),locLxSite(0))=1.0;
    // rdf2(parity(0),spinRow(0),locEoSite(0))=1.0;
  }
  verbosity_lv=1;
  
//  crash("");
//  
//  auto u=mergedComp<CompsList<Spin,ComplId>>(0);
//  
//  /////////////////////////////////////////////////////////////////
//  
//  master_printf("allocated in %p\n",conf._data);
//  {
//    auto c=conf.getWritable();
//    master_printf("allocated in %p\n",c._data);
//    double& e=c[locVol-1][3][2][2][1];
//    master_printf("end: %p, should be %p\n",&e,c._data+locVol*4*3*3*2-1);
//  }
//  
//  PAR(0,locVol,
//      CAPTURE(TO_WRITE(conf)),
//      ivol,
//      {
//	su3_put_to_id(conf[ivol][0]);
//      });
  
  // start_loc_rnd_gen(235235);
  
  // spincolor *in=nissa_malloc("in",locVol+bord_vol,spincolor);
  // spincolor *out=nissa_malloc("out",locVol+bord_vol,spincolor);
  // spincolor *tmp=nissa_malloc("tmp",locVol+bord_vol,spincolor);
  
  // /// First test: load a spincolor and unload it
  // generate_undiluted_source(in,RND_Z2,-1);
  
  // quda_iface::remap_nissa_to_quda(tmp,in);
  // quda_iface::remap_quda_to_nissa(out,tmp);
  
  // master_printf("testing map and unmap, residue: %lg\n",rel_diff_norm(out,in));
  
  // /// Second test: apply the dirac operator
  
  // quad_su3 *conf=nissa_malloc("conf",locVol+bord_vol+edge_vol,quad_su3);
  // spincolor *out_nissa=nissa_malloc("out_nissa",locVol+bord_vol,spincolor);
  
  // generate_hot_lx_conf(conf);
  // master_printf("plaq: %lg\n",global_plaquette_lx_conf(conf));
  
  // vector_reset(in);
  // in[0][0][0][0]=1.0;
  
  // const double kappa=0.1325,csw=1.345,mu=0.243;
  // quda_iface::apply_tmD(out,conf,kappa,csw,mu,in);
  
  // clover_term_t *Cl=nissa_malloc("Cl",locVol,clover_term_t);
  // inv_clover_term_t *invCl=nissa_malloc("invCl",locVol,inv_clover_term_t);
  // clover_term(Cl,csw,conf);
  // invert_twisted_clover_term(invCl,mu*tau3[0],kappa,Cl);
  // apply_tmclovQ(out_nissa,conf,kappa,Cl,mu,in);
  // nissa_free(invCl);
  // nissa_free(Cl);
  
  // safe_dirac_prod_spincolor(out_nissa,base_gamma[5],out_nissa);
  
  // master_printf("comparing\n");
  // for(int ivol=0;ivol<locVol;ivol++)
  //   for(int id=0;id<NDIRAC;id++)
  //     for(int ic=0;ic<NCOL;ic++)
  // 	for(int ri=0;ri<2;ri++)
  // 	  {
  // 	    const double n=out_nissa[ivol][id][ic][ri];
  // 	    const double q=out[ivol][id][ic][ri];
  // 	    if(fabs(n)>1e-10 or fabs(q)>1e-10)
  // 	      master_printf("out,[nissa,quda][ivol=%d,id=%d,ic=%d,ri=%d]: %lg %lg\n",
  // 			    ivol,id,ic,ri,n,q);
  // 	  }
  // master_printf("testing tmD, residue: %lg\n",rel_diff_norm(out,out_nissa));
  
  // nissa_free(tmp);
  // nissa_free(in);
  // nissa_free(out);
  // nissa_free(out_nissa);
  // nissa_free(conf);
}

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  in_main(narg,arg);
  close_nissa();
  
  return 0;
}
