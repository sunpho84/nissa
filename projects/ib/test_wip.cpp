#include <nissa.hpp>

// #include "base/quda_bridge.hpp"

using namespace nissa;
StackTens<OfComps<ComplId>,double> fd;

// double rel_diff_norm(spincolor *test,spincolor *ref)
// {
//   double_vector_subtassign((double*)test,(double*)ref,locVol*sizeof(spincolor)/sizeof(double));
//   const double norm2_diff=double_vector_glb_norm2(test,locVol);
//   const double norm2_ref=double_vector_glb_norm2(ref,locVol);
//   const double res=sqrt(norm2_diff/norm2_ref);
  
//   return res;
// }

namespace nissa
{
  DECLARE_TRANSPOSABLE_COMP(Spin,int,NDIRAC,spin);
  DECLARE_TRANSPOSABLE_COMP(Color,int,NCOL,color);
  DECLARE_TRANSPOSABLE_COMP(Fuf,int,1,fuf);
}

template <typename F,
	  typename Float=typename F::Fund>
Float plaquette(const F& conf)
{
  ASM_BOOKMARK_BEGIN("plaq");
  
  Field2<OfComps<ColorRow,ColorCln,ComplId>,Float> temp1,temp2;
  Field2<OfComps<>,Float> squares(0.0);
  
  for(Dir mu=0;mu<NDIM;mu++)
    for(Dir nu=mu+1;nu<NDIM;nu++)
      {
	temp1=(conf(mu)*shift(conf(nu),bw,mu));
	temp2=(conf(nu)*shift(conf(mu),bw,nu));
	
	squares+=real(trace(temp1*dag(temp2)));
      }
  
  const Float plaq=
    squares.glbReduce()/glbVol/(2*NCOL*NCOL);
  
  ASM_BOOKMARK_END("plaq");
  
  return plaq;
}

namespace nissa
{
  DECLARE_UNTRANSPOSABLE_COMP(FFTComp,int64_t,0,fftComp);
  DECLARE_UNTRANSPOSABLE_COMP(FFTOrthoComp,int64_t,0,fftOtherComp);
  //DECLARE_UNTRANSPOSABLE_COMP(FFTInternalDegComp,int64_t,0,fftinternalDegComp);
  
  template <typename T>
  struct FFT
  {
    using InternalDegComps=
      TupleFilterAllTypes<typename T::Comps,
			  CompsList<ComplId,LocLxSite>>;
    
    using FFTInternalDegComp=
      MergedComp<InternalDegComps>;
    
    using FFTMergedComp=
      MergedComp<CompsList<FFTOrthoComp,FFTInternalDegComp>>;
    
    using BufComps=
      OfComps<FFTComp,FFTMergedComp>;
    
    using BufEl=
      MergedComp<BufComps>;
    
    using Fund=
      typename T::Fund;
    
    using Compl=
      std::array<Fund,2>;
    
    struct Sizes
    {
      const int dir;
      
      const FFTComp dirSize;
      
      const FFTOrthoComp orthoDirSize;
      
      const FFTInternalDegComp nInternalDegs;
      
      const FFTMergedComp mergedCompSize;
      
      Sizes(const T& t,
	    const int dir) :
	dir{dir},
	dirSize{locSize[dir]},
	orthoDirSize{t.nSites()()/locSize[dir]},
	nInternalDegs{t.nInternalDegs/2},
	mergedCompSize{nInternalDegs()*orthoDirSize()}
      {
      }
      
      INLINE_FUNCTION CUDA_HOST_AND_DEVICE
      std::tuple<FFTComp,FFTOrthoComp> decomposeCoords(const coords_t& coords) const
      {
	const FFTComp fftComp=
	  coords[dir];
	
	FFTOrthoComp fftOrthoComp=0;
	
	for(int iDir=0;iDir<NDIM-1;iDir++)
	  {
	    const Dir orthoDir=perp_dir[dir][iDir];
	    fftOrthoComp=fftOrthoComp*glbSize[orthoDir()]+coords[orthoDir()];
	  }
	
	return {fftComp,fftOrthoComp};
      }
      
      INLINE_FUNCTION CUDA_HOST_AND_DEVICE
      std::tuple<FFTComp,FFTOrthoComp> decomposeSite(const int& site) const
      {
	return decomposeCoords(locCoordOfLoclx[site]);
      }
      
      INLINE_FUNCTION CUDA_HOST_AND_DEVICE
      std::tuple<coords_t,FFTInternalDegComp> getCoordsInternalDegFromBufComps(const FFTComp& fftComp,
									       const FFTMergedComp& mergedComp) const
      {
	auto [fftOrthoComp,fftInternalDegComp]=
	  mergedComp.decompose(std::make_tuple(orthoDirSize));
	
	coords_t coords;
	coords[dir]=fftComp();
	
	for(int iDir=NDIM-1;iDir>=0;iDir--)
	  {
	    const Dir orthoDir=perp_dir[dir][iDir];
	    coords[orthoDir()]=fftOrthoComp()%glbSize[orthoDir()];
	    fftOrthoComp=fftOrthoComp/glbSize[orthoDir()];
	  }
	
	return {coords,fftInternalDegComp};
      }
    };
    
    static void fft(T& f)
    {
      const int dir0=0;
      
      Sizes sizes(f,dir0);
      
      DynamicTens<BufComps,Compl,T::execSpace> buf(std::make_tuple(sizes.dirSize,sizes.mergedCompSize));
      
      PAR(0,locVol,CAPTURE(sizes,
			   TO_WRITE(buf),
			   TO_READ(f)),
	  site,
	  {
	    compsLoop<InternalDegComps>([&sizes,
					 &buf,
					 &site,
					 &f](const auto&...c)
	    {
	      const auto [fftComp,fftOrthoComp]=
		sizes.decomposeSite(site);
	      
	      const FFTMergedComp mergedComp(std::make_tuple(sizes.orthoDirSize),fftOrthoComp,c...);
	      
	      for(int rI=0;rI<2;rI++)
		buf(fftComp,mergedComp)[rI]=f(LocLxSite(site),c...,reIm(rI));
	    },std::make_tuple());
	  });
      
      /////////////////////////////////////////////////////////////////
      
      for(int dirFrom=0;dirFrom<NDIM-1;dirFrom++)
	{
	  const int dirTo{dirFrom+1};
	  
	  const Sizes sizesFrom(f,dirFrom);
	  
	  const Sizes sizesTo(f,dirTo);
	  
	  DynamicTens<BufComps,Compl,T::execSpace> buf2(std::make_tuple(sizesTo.dirSize,sizesTo.mergedCompSize));
	  
	  PAR(0,sizesFrom.mergedCompSize()*sizesFrom.dirSize(),
	      CAPTURE(sizesFrom,
		      sizesTo,
		      TO_WRITE(buf2),
		      TO_READ(buf)),
	      elemFrom,
	      {
		const auto [fftCompFrom,mergedCompFrom]=
		  BufEl::decompose(std::make_tuple(sizesFrom.dirSize,sizesFrom.mergedCompSize),elemFrom);
		
		const auto [siteCoords,internalDeg]=
		  sizesFrom.getCoordsInternalDegFromBufComps(fftCompFrom,mergedCompFrom);
		
		const auto [fftCompTo,fftOrthoCompTo]=
		  sizesTo.decomposeCoords(siteCoords);
		
		const FFTMergedComp mergedCompTo(std::make_tuple(sizesTo.orthoDirSize),fftOrthoCompTo);
		
		for(int rI=0;rI<2;rI++)
		  buf2(fftCompTo,mergedCompTo)[rI]=
		    buf(fftCompTo,FFTMergedComp(mergedCompFrom))[rI];
	      });
	  
	  std::swap(buf,buf2);
	}
      
      /////////////////////////////////////////////////////////////////
      
      // Sizes sizes3(f,3);
      // PAR(0,sizes3.mergedCompSize(),
      // 	  CAPTURE(sizes3,
      // 		  TO_READ(buf),
      // 		  TO_WRITE(f)),
      // 	  mergedComp,
      // 	  {
      // 	    // const auto [siteCoords,internalDeg]=
      // 	    //   sizes3.decomposeMergedComp(mergedComp);
	    
      // 	  //   for(FFTComp fftComp=0;fftComp<glbSize[3];fftComp++)
      // 	  //     {
      // 	  // 	coords_t c;
      // 	  // 	for(int iDir=0;iDir<NDIM-1;iDir++)
      // 	  // 	  {
      // 	  // 	    const Dir orthoDir=perp_dir[3][iDir];
      // 	  // 	    fftOrthoComp=fftOrthoComp*glbSize[orthoDir()]+coords[orthoDir()];
      // 	  // }
      // 	  // 	for()
      // 	  //   internalDeg.decompose();
      // 	    //   const FFTMergedComp mergedComp=
      // 	    // 	index(std::make_tuple(sizes.orthoDirSize),fftOrthoComp,c...);
	      
      // 	    //   for(int rI=0;rI<2;rI++)
      // 	    // 	buf(fftComp,mergedComp)[rI]=f(LocLxSite(site),c...,reIm(rI));
      // 	    // },std::make_tuple());
      // 	  });
      
    }
  };
}

void in_main(int narg,char **arg)
{
  const int T=8,L=4;
  
  init_grid(T,L);
  
  //constexpr Float128 a=(1e-60);
  //const double c=(a<0.5)?std::log1p(-a.roundDown()):log((1-a).roundDown());
  
  // for(Float128 b=(Float128)0+(1e-60);b<=1;b*=1.1)
  //   {
  //     const Float128 c=(0.25+b)*2;
  //     const double d=
  // 	((c>=1.0/4 and c<3.0/4) or (c>=5.0/4 and c<7.0/4))?
  // 	sin(M_PI*(0.5-c).roundUp()):
  // 	cos(M_PI*c.roundDown());
  //     printf("%.017lg %.017lg  %.017lg %.017lg\n",b.roundDown(),b.roundDown()*2*M_PI,d,cos(M_PI*c.roundDown()));
  //   }

  // auto uu=[](const uint64_t& seed)
  // {
  //   Rng r(seed);
  //   auto u=[&r](const int i)
  //   {
  //     return r.encrypter.key[i];
  //   };
  //   printf("%lu %lu %lu %lu %lu\n",u(0),u(1),u(2),u(3),u(4));
  // };

  // uu(0);
  // uu(1);
  // uu(353252626);
  // uu(353252627);
  
  // StackTens<OfComps<Fuf>,double> er;
  // er=0.0;
  // const double ee=er;
  
  RngState state(32534643);
  printf("%lu\n",state.counter[0]);
  const auto view=state.getNormalDistr(2);
  printf("%lu\n",state.counter[0]);
  const auto view2=state.getNormalDistr(2);
  printf("%lu\n",state.counter[0]);
  PAR(0,1,CAPTURE(view,view2),i,
      {
	printf("%lg\n",view.draw(i));
	printf("%lg\n",view2.draw(i));
      });

  constexpr uint32_t bbb=(1u<<31)+(1u<<30);
  const double d=ProbDistr::NormalRngDistr::transformCos({~0u-1,bbb-1});
  const double e=ProbDistr::NormalRngDistr::transformCos({~0u,bbb-1});
  const double f=ProbDistr::NormalRngDistr::transformCos({0u,bbb});
  const double g=ProbDistr::NormalRngDistr::transformCos({1u,bbb});
  printf("d d d %.017lg\n",d);
  printf("d d d %.017lg\n",e);
  printf("d d d %.017lg\n",f);
  printf("d d d %.017lg\n",g);
  
  // printf("%lu\n",state.counter.val[0]);
  // RngGaussDistrView gaussDistr(state,2);
  // printf("%lu\n",state.counter.val[0]);
  // printf("%lg\n",gaussDistr.draw(0));
  // printf("%lg\n",gaussDistr.draw(1));
  // printf("%lg\n",gaussDistr.draw(2));
  // double min=1;
  // for(uint64_t i=0;i<100000000000;i++)
  //   {
  //     double x=Rng::transform32bitsListIntoUniform(i,i>>32).roundDown();
  //     // if(x<min)
  //     // 	{
  //     printf("%lu sss %lg\n",i,x);
  // 	//   min=x;
  // 	// }
  //   }
  
  //return ;
  LxField<quad_su3> conf("conf");
  {
    DynamicTens<OfComps<SpinRow,LocEoSite>,double, MemoryType::CPU> d(std::make_tuple(locEoSite(locVolh)));
    auto e=d.getReadable();
    // e(spinRow(1))=1.0;
  }
  
  {
    // verbosity_lv=3;
    // Field2<OfComps<SpinRow,ComplId>,double> df(WITH_HALO);
    
    // for(int loclxsite=0;loclxsite<locVol;loclxsite++)
    //   df(LocLxSite(loclxsite),Spin(0),Re)=loclxsite;
    // printf("%lg %d\n",shift(conj(df),bw,Dir(2))(LocLxSite(0),Spin(0),Re),loclxNeighup[0][2]);
    // printf("%lg %d\n",shift(df,fw,Dir(2))(LocLxSite(0),Spin(0),Re),loclxNeighdw[0][2]);
    
    LxField<quad_su3> gcr("gcr",WITH_HALO);
    read_ildg_gauge_conf(gcr,"../test/data/L4T8conf");
    Field2<OfComps<Dir,ColorRow,ColorCln,ComplId>,double,FULL_SPACE,FieldLayout::GPU,MemoryType::CPU> gcH(WITH_HALO);
    
    gcH.locLxSite(2);
    ASM_BOOKMARK_BEGIN("expI");
    StackTens<OfComps<Dir>,double> phases{1,0,0,0};
    StackTens<OfComps<Dir>,int> glbSizes{8,4,4,4};
    StackTens<OfComps<Dir,ComplId>,double> expI=exp(I*M_PI*phases/glbSizes);
    ASM_BOOKMARK_END("expI");
    
    master_printf("Let's see if we complain: %lg %lg\n",real(I),imag(I));
    master_printf("Let's see if we complain: %lg %lg\n",real(expI)(Dir(0)),imag(expI)(Dir(0)));
    
    master_printf("Copying the conf to field2\n");
    for(LocLxSite site=0;site<locVol;site++)
      for(Dir mu=0;mu<NDIM;mu++)
	for(ColorRow cr=0;cr<NCOL;cr++)
	  for(ColorCln cc=0;cc<NCOL;cc++)
	    for(ComplId reIm=0;reIm<2;reIm++)
	      gcH(site,mu,cr,cc,reIm)=gcr[site()][mu()][cr()][cc()][reIm()];
    Field2<OfComps<Dir,ColorRow,ColorCln,ComplId>,double> gc(WITH_HALO);
    master_printf("Copying the conf from host to device\n");
    gc=gcH;
    master_printf("plaq with new method: %.16lg\n",plaquette(gc));
    Field2<OfComps<Dir,ColorRow,ColorCln,ComplId>,double> gcp(WITH_HALO);
    
    // auto rrr=(gc*expI).template closeAs<std::decay_t<decltype(gc.createCopy())>>();
    gc*=expI;
    
    auto e=(gc*expI).closeAs(gc);
    
    master_printf("plaq after phasing: %.16lg\n",plaquette(gc));
    master_printf("plaq after phasing: %.16lg\n",plaquette(e));
    FFT<std::decay_t<decltype(e)>>::fft(e);
    
// memcpy(gc.data.storage,gcr._data,sizeof(quad_su3)*locVol);
    // Field2<OfComps<>,double> t(WITH_HALO);
    // Field2<OfComps<>,double> u;
    // t=0.875;
    // u=shift(t,fw,Dir(0))*shift(t,fw,Dir(0));
    // printf("u: %lg\n",u(LocLxSite(0)));
    
    
    // df.updateHalo();
    // auto fl=df.flatten();
    // fl(decltype(df)::FlattenedInnerComp(0));
    // Field2<OfComps<ComplId>,double> ds;
    // // EoField2<OfComps<SpinRow>,double> df2;
    // //auto rdf=df.getWritable();
    // // auto rdf2=df2.getWritable();
    // df=I;
    // master_printf("written 0? %lg\n",df(locLxSite(0),spinRow(0),reIm(0)));
    // master_printf("written 1? %lg\n",df(locLxSite(0),spinRow(0),reIm(1)));
    // df=df+I;
    // ds=trace(dag(df(spinRow(0))));
    // master_printf("written 2? %lg\n",df(locLxSite(0),spinRow(0),reIm(1)));
    // master_printf("copied -2? %lg\n",ds(locLxSite(0),reIm(1)));
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
