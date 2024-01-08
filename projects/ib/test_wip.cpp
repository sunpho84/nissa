#include <memory>

#include <nissa.hpp>

#include <tuples/tuple.hpp>

using namespace nissa;


void ddd()
{
  static_assert(std::is_move_constructible_v<const CompsBinder<std::tuple<LocLxSite>,
		const MirroredNode<std::tuple<LocLxSite, Dir>, DynamicTens<std::tuple<LocLxSite, Dir>, const GlbCoord, MemoryType::CPU, true>, DynamicTens<std::tuple<LocLxSite, Dir>, const GlbCoord, MemoryType::GPU, true>, const GlbCoord> , std::tuple<Dir>, const GlbCoord>>,"ccccc");
  }

template <typename F,
	  typename Float=typename F::Fund>
Float plaquette(const F& conf)
{
  ASM_BOOKMARK_BEGIN("plaq");
  
  Field<OfComps<ColorRow,ColorCln,ComplId>,Float> temp1,temp2;
  Field<OfComps<>,Float> squares(0.0);
  
  for(Dir mu=0;mu<NDIM;mu++)
    for(Dir nu=mu+1;nu<NDIM;nu++)
      {
	temp1=(conf(mu)*shift(conf(nu),bw,mu));
	temp2=(conf(nu)*shift(conf(mu),bw,nu));
	
	squares+=real(trace(temp1*dag(temp2)));
      }
  
  const Float plaq=
    squares.glbReduce()/lat->getGlbVol()/(2*NCOL*NCOL);
  
  ASM_BOOKMARK_END("plaq");
  
  return plaq;
}

namespace nissa
{
  /// One of the fourth roots of unity
  enum class IdFourthRoot{One,I,MinusOne,MinusI};
  
  /// Product of two fourth roots of unity
  constexpr IdFourthRoot operator*(const IdFourthRoot& first,
				   const IdFourthRoot& second)
  {
    return IdFourthRoot(((int)first+(int)second)%4);
  }
  
  /// Complex conjugate of one of the fourth roots of unity
  constexpr IdFourthRoot operator~(const IdFourthRoot& i)
  {
    return IdFourthRoot(((int)i+((int)i%2)*2)%4);
  }
  
  /// Negation of one of the fourth roots of unity
  constexpr IdFourthRoot operator-(const IdFourthRoot& i)
  {
    return IdFourthRoot(((int)i+2)%4);
  }
  
  struct Gamma
  {
    using Entry=
      std::pair<SpinCln,IdFourthRoot>;
    
    using Entries=
      StackTens<CompsList<SpinRow>,Entry>;
    
    const Entries entries;
    
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    constexpr Gamma operator*(const Gamma& oth) const
    {
      Entries tmp;
      
      for(SpinRow ig1=0;ig1<nDirac;ig1++)
	{
	  //This is the line to be taken on the second matrix
	  const SpinRow ig2=transp(entries(ig1).first);
	  
	  //The entries of the output is, on each line, the complex
	  //product of the entries of the first matrix on that line, for
	  //the entries of the second matrix on the line with the index
	  //equal to the column of the first matrix which is different
	  //from 0 (which again is ig2)
	  const IdFourthRoot r=
	    entries(ig1).second*oth.entries(ig2).second;
	  
	  //For each line, the column of the output matrix which is
	  //different from 0 is the column of the second matrix different
	  //from 0 on the line with index equal to the column of the first
	  //matrix which is different from 0 (that is, ig2)
	  tmp(ig1)={oth.entries(ig2).first,r};
	};
      
      return {tmp};
    }
    
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    constexpr Gamma operator-() const
    {
      Entries tmp;
      
      for(SpinRow i=0;i<nDirac;i++)
	tmp(i)={tmp(i).first,-tmp(i).second};
      
      return {tmp};
    }
    
    /// Returns the hermitian
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    constexpr Gamma dag() const
    {
      Entries tmp;
      
      for(SpinRow sr=0;sr<nDirac;sr++)
	{
	  const auto& [sc,c]=this->entries(sr);
	  tmp(transp(sc))={sc,~c};
	}
      
      return {tmp};
    }
  };
  
  //static constexpr Gamma Id{{{0,IdFourthRoot::One}}};
}

namespace GammaBasis
{
#define E(P,V) Gamma::Entry{P,IdFourthRoot::V}
#define GAMMA(NAME,P0,V0,P1,V1,P2,V2,P3,V3)				\
  constexpr Gamma NAME{{E(P0,V0),E(P1,V1),E(P2,V2),E(P3,V3)}}
  
  GAMMA(X,3,MinusI,2,MinusI,1,I,0,I);
  GAMMA(Y,3,MinusOne,2,One,1,One,0,MinusOne);
  GAMMA(Z,2,MinusI,3,I,0,I,1,MinusI);
  GAMMA(T,2,MinusOne,3,MinusOne,0,MinusOne,1,MinusOne);
  
#undef GAMMA
#undef E
  constexpr Gamma G5=T*X*Y*Z;
  constexpr Gamma e=-G5.dag();
  constexpr std::array<Gamma,16> basis{X*X,X,Y,Z,T,G5,G5*X,G5*Y,G5*Z,G5*T,T*X,T*Y,T*Z,Y*Z,Z*X,X*Y};
}


using PhotonField=
  Field<CompsList<Dir>>;

using PhotonProp=
  Field<CompsList<DirRow,DirCln>>;

using ComplVectorProp=
  Field<CompsList<DirRow,DirCln,ComplId>>;

using ComplVectorField=
  Field<CompsList<Dir,ComplId>>;

using ScalarField=
  Field<CompsList<>>;

PhotonField getZ2PhotonField(const int offset=0)
{
  PhotonField res;
  
  RngState rng(235235+offset);
  
  // auto printRng=[](const RngState& r)
  // {
  //   for(const auto& i : r.counter)
  //     master_printf("%u\n",i);
  // };
  
  // printRng(rng);
  // double e=take_time();
  // for(int i=0;i<300;i++)
  //   {
  auto d=rng.getZ2Distr(lat->getGlbVol()()*nDim);
  
  PAR_ON_EXEC_SPACE(ScalarField::execSpace,
		    0,
		    lat->getLocVol(),
		    CAPTURE(TO_WRITE(res),
			    c=lat->getGlbLxOfLocLx().getRef(),
			    d),
		    site,
			{
			  ASM_BOOKMARK_BEGIN("rnd");
			  for(Dir dir=0;dir<nDim;dir++)
			    res(site,dir)=d.draw(dir()+nDim()*c(site)());
			  ASM_BOOKMARK_END("rnd");
			});
  // printRng(rng);
      
      // master_printf("ciao %lg\n",res.copyToMemorySpaceIfNeeded<MemoryType::CPU>()(LocLxSite(0)));
  //   }
  // master_printf("ciao took %lg s\n",take_time()-e);
  
  return res;
}

void in_main(int narg,char **arg)
{
  
  // constexpr ExecSpace ex(MemoryType::CPU);
  // constexpr MemoryType mt=getMemoryType<ex>();
  // constexpr bool t=ex.template runOn<MemoryType::CPU>();
  // constexpr int m=(execOnCPU*execOnGPU).hasUniqueExecSpace();
  const int T=8,L=4;
  
  _lat=new LatticeResources;
  _lat->init(T,L);
  lat=std::make_unique<Lattice>(_lat->getRef());
  
#ifdef USE_QUDA
    if(use_quda) quda_iface::initialize();
#endif
     
  localizer::init();
  initFftw();
  
  /////////////////////////////////////////////////////////////////

   masterPrintf("TEST\n");
    // cudaGenericKernel<<<gridDimension,blockDimension>>>(0,locVol,[P=P.getWritable(),l=lat->glbCoordsOfLocLx.getReadable()] CUDA_DEVICE(auto i) mutable{P(LocLxSite(i))=l(LocLxSite(i));});
    // 	decrypt_cuda_error(cudaDeviceSynchronize(),"during kernel executionssss");
  
   PhotonField P=2*sin(M_PI*(lat->getGlbCoordsOfLocLx()+std::numeric_limits<double>::epsilon())/lat->getGlbSizes());
  // auto o=P.copyToMemorySpaceIfNeeded<MemoryType::CPU>().locLxSite(23).dirRow(3);
  // printf("AAA %ld %lg\n",ii(),o);
  
#define NISSA_HOST_DEVICE_LAMBDA(CAPTURES,BODY...)	\
  std::make_tuple([CAPTURE(CAPTURES)] BODY,[CAPTURE(CAPTURES)] DEVICE_ATTRIB BODY)
   
  P.onEachSite(NISSA_HOST_DEVICE_LAMBDA(CAPTURE(),
					(auto& P,const LocLxSite& site)
					{
					  P(site);
					}));
  
  const auto Ptilde=
    P*Lattice::perpDirs[Dir(0)];
  
   const ScalarField P2=
     dag(P)*P;

   const ScalarField zmt=(1-lat->spatialOriginsMask());
   printf("AAA %lg\n",zmt.copyToMemorySpaceIfNeeded<MemoryType::CPU>().locLxSite(0));
   
   const ScalarField P2tilde=
     dag(Ptilde)*Ptilde;
   
   PhotonProp test=kronDelta<DirRow,DirCln>-(P*dag(Ptilde)+Ptilde*dag(P)-P*dag(P))/P2tilde;
   //test.dirRow(0).dirCln(0)=0;
   
   Field<OfComps<Dir,ColorRow,ColorCln,ComplId>,double> tmp;

   readFieldFromIldgFile(tmp,"/tmp/conf","ildg-binary-data");
   
   for(LocLxSite site=0;site<lat->getLocVol();site++)
     {
       for(Dir dir=0;dir<4;dir++)
	 {
	   double f=tmp.locLxSite(site).dirRow(dir).colorCln(0).colorRow(0).reIm(0);
	   double g=lat->getGlbCoordsOfLocLx(site,dir)();
	   if(f!=g)
	     CRASH("not agreeing");
	 }
     }
   readFieldFromIldgFile(tmp,"../test/data/L4T8conf","ildg-binary-data");
   masterPrintf("plaq with new method: %.16lg\n",plaquette(tmp));

   
   CRASH("");
   auto printPhotonProp=
     [](const PhotonProp& prop)
     {
       for(DirRow r=0;r<nDim;r++)
	 {
	   // for(DirCln c=0;c<nDim();c++)
	   //   printf("%lg ",prop.locLxSite(glblx_of_coord({3,1,2,1}))(r,c));
	   // printf("\n");
	 }
       printf("\n");
     };
   
   printPhotonProp(test);
   
   test.dirRow(0).dirCln(0)=sqrt(test.dirRow(0).dirCln(0));
   const PhotonProp test2=test*test;
   printPhotonProp(test2);
   
   const PhotonProp prop=
     (1-lat->spatialOriginsMask())*
     (kronDelta<DirRow,DirCln>-(P*dag(Ptilde)+Ptilde*dag(P)-P*dag(P))/P2tilde)/P2;
  
  const auto eta=
    getZ2PhotonField();
  
  const PhotonField phi=
    real(fft<ComplVectorField>(-1,prop*fft<ComplVectorField>(+1,eta)));
  
  printf("AAA %lg\n",phi.copyToMemorySpaceIfNeeded<MemoryType::CPU>().locLxSite(23).dirRow(0));
  CRASH("");
  /////////////////////////////////////////////////////////////////
  
  
  Field<CompsList<ColorRow>> r;
  Field<CompsList<ColorRow,ComplId>> rc;
  fft(rc,+1,r);

    Field<CompsList<ColorRow,ComplId>> e;
    e=conj(e);
  
  // const int s=glblx_of_coord({1,0,0,0});
  // PAR(0,
  //     lat->getLocVol(),
  //     CAPTURE(TO_WRITE(e),
  // 	      s,
  // 	      g=lat->getGlbLxOfLocLx().getRef()),
  //     site,
  //     {
  // 	for(ColorRow cr=0;cr<3;cr++)
  // 	  for(ComplId ri=0;ri<2;ri++)
  // 	    e(site,cr,ri)=g[site]==s and ri==0;
  //     });

  // decltype(auto) fDone=
  //   fft<Field<CompsList<ColorRow,ComplId>>>(+1,e).copyToMemorySpaceIfNeeded<MemoryType::CPU>();
    
  CRASH("uncomment");
  // for(LocLxSite site=0;
  //     site<lat->getLocVol();
  //     site++)
  //   {
  //     auto p=fDone.locLxSite(site);
  //     auto c=p.colorRow(0);
  //     auto g=glb_coord_of_glblx(lat->getGlbLxOfLocLx(site));
  //     master_printf("site %ld (%d,%d,%d,%d): %lg %lg\n",
  // 		    site(),g[0],g[1],g[2],g[3],c.reIm(0),c.reIm(1));

  //     // double r=0;
  //     // for(Dir d=0;d<nDim;d++)
  //     // 	r+=M_PI*glbCoordOfLoclx[site()][d()];
      
  //     for(ColorRow cr=0;cr<3;cr++)
  // 	for(ComplId ri=0;ri<2;ri++)
  // 	  if(p(cr,ri)!=c(ri))
  // 	    master_printf("site %ld cr %d ri %d: %lg %lg vs %lg %lg\n",site(),cr(),ri(),p(cr).reIm(0),p(cr).reIm(1),c.reIm(0),c.reIm(1));
	    
  //   }
  
  

  // DynamicTens<OfComps<SpinRow,LocEoSite>,double,MemoryType::CPU>
  // din(std::make_tuple(locEoSite(locVolh))); din=0;

  //     using MC=MergedComp<CompsList<SpinRow,LocEoSite>>;
  //     auto
  //     mc=MC::merge(std::make_tuple(LocEoSite(locVolh)),spinRow(2),locEoSite(34));
  //     ein=124;
  //     ein(mc)=23124;
  //     master_printf("%lg\n",din.spinRow(2).locEoSite(34));

  //     invokeWithTypesOfTuple<CompsList<SpinRow,LocEoSite>>([&din]<typename...T>()
  // 							 {
  // 							   ((master_printf("original comp %s size
  // %d\n",demangle(typeid(T).name()).c_str(),din.getCompSize<T>())),...);
  // 							 });
  //     invokeWithTypesOfTuple<decltype(ein)::DynamicComps>([&ein]<typename
  //     T>()
  // 						      {
  // 							master_printf("merged comp %s size
  // %d\n",demangle(typeid(T).name()).c_str(),ein.getCompSize<T>());
  // 						      });

  //     {
  AllToAllComm<GlbLxSite, LocLxSite> aa;

  auto ff = [](const LocLxSite &locLxSite) -> CompsList<MpiRank, GlbLxSite> {
    return {0, lat->getGlbLxOfLocLx(locLxSite)()};
  };
  masterPrintf("%s\n", demangle(typeid(ff).name()).c_str());

  constexpr MemoryType execSpace = MemoryType::
#ifdef USE_CUDA
      GPU
#else
      CPU
#endif
      ;
      
  DynamicTens<OfComps<SpinRow,LocLxSite,ColorCln>,double,MemoryType::CPU> _in(std::make_tuple(lat->getLocVol()));
       
      aa.init(lat->getLocVol(),
 	      ff);
      
      for(LocLxSite locLxSite=0;locLxSite<lat->getLocVol();locLxSite++)
	_in(locLxSite)=lat->getGlbLxOfLocLx(locLxSite)();
      
      decltype(auto) in=_in.template copyToMemorySpaceIfNeeded<execSpace>();
      auto _out=aa.communicate(in);
      
      decltype(auto) out=_out.template copyToMemorySpaceIfNeeded<MemoryType::CPU>();
      if(isMasterRank())
	for(GlbLxSite glbLxSite=0;glbLxSite<lat->getGlbVol();glbLxSite++)
	  if(const double o=out(glbLxSite).spinRow(0).colorCln(0),g=glbLxSite();o!=g)
	    CRASH("%lg %lg\n",o,g);
      masterPrintf("alright!\n");
      
      auto bb=aa.inverse();
      DynamicTens<OfComps<SpinRow,LocLxSite,ColorCln>,double,execSpace> _back(std::make_tuple(lat->getLocVol()));
      bb.communicate(_back,_out);
      decltype(auto) back=_back.template copyToMemorySpaceIfNeeded<MemoryType::CPU>();
      
      for(LocLxSite locLxSite=0;locLxSite<lat->getLocVol();locLxSite++)
	if(const double b=back(locLxSite).spinRow(0).colorCln(0),g=lat->getGlbLxOfLocLx(locLxSite)();b!=g)
	    CRASH("%lg %lg\n",b,g);
      masterPrintf("alright2!\n");
      
      auto cc=bb*aa;
      decltype(auto) shouldBeId=cc.communicate(in).template copyToMemorySpaceIfNeeded<MemoryType::CPU>();
      for(LocLxSite locLxSite=0;locLxSite<lat->getLocVol();locLxSite++)
	if(const double s=shouldBeId(locLxSite).spinRow(0).colorCln(0),
	   i=_in(locLxSite).spinRow(0).colorCln(0);
	   s!=i)
	    CRASH("%lg %lg\n",s,i);
      masterPrintf("alright3!\n");
      
      
//     //int aa=(time(0)%2)?[](){printf("This\n");return 0;}():[](){printf("That\n");return 1;}();
    
//   return ;
//   //constexpr Float128 a=(1e-60);
//   //const double c=(a<0.5)?std::log1p(-a.roundDown()):log((1-a).roundDown());
  
//   // for(Float128 b=(Float128)0+(1e-60);b<=1;b*=1.1)
//   //   {
//   //     const Float128 c=(0.25+b)*2;
//   //     const double d=
//   // 	((c>=1.0/4 and c<3.0/4) or (c>=5.0/4 and c<7.0/4))?
//   // 	sin(M_PI*(0.5-c).roundUp()):
//   // 	cos(M_PI*c.roundDown());
//   //     printf("%.017lg %.017lg  %.017lg %.017lg\n",b.roundDown(),b.roundDown()*2*M_PI,d,cos(M_PI*c.roundDown()));
//   //   }

//   // auto uu=[](const uint64_t& seed)
//   // {
//   //   Rng r(seed);
//   //   auto u=[&r](const int i)
//   //   {
//   //     return r.encrypter.key[i];
//   //   };
//   //   printf("%lu %lu %lu %lu %lu\n",u(0),u(1),u(2),u(3),u(4));
//   // };

//   // uu(0);
//   // uu(1);
//   // uu(353252626);
//   // uu(353252627);
  
//   // StackTens<OfComps<Fuf>,double> er;
//   // er=0.0;
//   // const double ee=er;
  
//   RngState state(32534643);
//   printf("%lu\n",state.counter[0]);
//   const auto view=state.getNormalDistr(2);
//   printf("%lu\n",state.counter[0]);
//   const auto view2=state.getNormalDistr(2);
//   printf("%lu\n",state.counter[0]);
//   PAR(0,1,CAPTURE(view,view2),i,
//       {
// 	printf("%lg\n",view.draw(i));
// 	printf("%lg\n",view2.draw(i));
//       });

//   constexpr uint32_t bbb=(1u<<31)+(1u<<30);
//   const double d=ProbDistr::NormalRngDistr::transformCos({~0u-1,bbb-1});
//   const double e=ProbDistr::NormalRngDistr::transformCos({~0u,bbb-1});
//   const double f=ProbDistr::NormalRngDistr::transformCos({0u,bbb});
//   const double g=ProbDistr::NormalRngDistr::transformCos({1u,bbb});
//   printf("d d d %.017lg\n",d);
//   printf("d d d %.017lg\n",e);
//   printf("d d d %.017lg\n",f);
//   printf("d d d %.017lg\n",g);

//   master_printf("/////////////////////////// defining //////////////////////////////////////\n");
//   MirroredNode<CompsList<LocLxSite,ComplId>,Field2<CompsList<ComplId>,double,FieldCoverage::FULL_SPACE,FieldLayout::CPU,MemoryType::CPU>> a(doNotAllocate);
//   master_printf("//////////////////////////// allocating /////////////////////////////////////\n");
//   a.allocate(WITH_HALO);
//   master_printf("//////////////////////////// assigning /////////////////////////////////////\n");
//   a=1;
//   master_printf("///////////////////////////// updating the copy ////////////////////////////////////\n");
//   a.updateDeviceCopy();
//   master_printf("/////////////////////////////////////////////////////////////////\n");
//   master_printf("Is 1? %lg\n",a.locLxSite(0).reIm(0));

//   // auto b=a.deviceVal.template copyToMemorySpace<MemoryType::CPU>();
//   // master_printf("Is 1? %lg\n",b.locLxSite(0).reIm(0));
//   // crash("");

//   // printf("%lu\n",state.counter.val[0]);
//   // RngGaussDistrView gaussDistr(state,2);
//   // printf("%lu\n",state.counter.val[0]);
//   // printf("%lg\n",gaussDistr.draw(0));
//   // printf("%lg\n",gaussDistr.draw(1));
//   // printf("%lg\n",gaussDistr.draw(2));
//   // double min=1;
//   // for(uint64_t i=0;i<100000000000;i++)
//   //   {
//   //     double x=Rng::transform32bitsListIntoUniform(i,i>>32).roundDown();
//   //     // if(x<min)
//   //     // 	{
//   //     printf("%lu sss %lg\n",i,x);
//   // 	//   min=x;
//   // 	// }
//   //   }

//   //return ;
//   LxField<quad_su3> conf("conf");
//   {
//     DynamicTens<OfComps<SpinRow,LocEoSite>,double,MemoryType::CPU> d(std::make_tuple(locEoSite(locVolh)));
//     auto dd=mergeComps<CompsList<SpinRow,LocEoSite>>(d);
//     auto e=d.mergeComps<CompsList<SpinRow,LocEoSite>>();
//     e=1;
    
//     // e(spinRow(1))=1.0;
//   }
  
//   {
//     // verbosity_lv=3;
//     // Field2<OfComps<SpinRow,ComplId>,double> df(WITH_HALO);
    
//     // for(int loclxsite=0;loclxsite<locVol;loclxsite++)
//     //   df(LocLxSite(loclxsite),Spin(0),Re)=loclxsite;
//     // printf("%lg %d\n",shift(conj(df),bw,Dir(2))(LocLxSite(0),Spin(0),Re),loclxNeighup[0][2]);
//     // printf("%lg %d\n",shift(df,fw,Dir(2))(LocLxSite(0),Spin(0),Re),loclxNeighdw[0][2]);
    
      //LxField<quad_su3> gcr("gcr",WITH_HALO);
//     Field2<OfComps<Dir,ColorRow,ColorCln,ComplId>,double,FULL_SPACE,FieldLayout::GPU,MemoryType::CPU> gcH(WITH_HALO);
    
//     gcH.locLxSite(2);
//     ASM_BOOKMARK_BEGIN("expI");
//     StackTens<OfComps<Dir>,double> phases{1,0,0,0};
//     StackTens<OfComps<Dir>,int> glbSizes{8,4,4,4};
//     StackTens<OfComps<Dir,ComplId>,double> expI=exp(I*M_PI*phases/glbSizes);
//     ASM_BOOKMARK_END("expI");
    
//     master_printf("Let's see if we complain: %lg %lg\n",real(I),imag(I));
//     master_printf("Let's see if we complain: %lg %lg\n",real(expI)(Dir(0)),imag(expI)(Dir(0)));
    
//     master_printf("Copying the conf to field2\n");
//     for(LocLxSite site=0;site<locVol;site++)
//       for(Dir mu=0;mu<NDIM;mu++)
// 	for(ColorRow cr=0;cr<NCOL;cr++)
// 	  for(ColorCln cc=0;cc<NCOL;cc++)
// 	    for(ComplId reIm=0;reIm<2;reIm++)
// 	      gcH(site,mu,cr,cc,reIm)=gcr[site()][mu()][cr()][cc()][reIm()];
    Field<OfComps<Dir,ColorRow,ColorCln,ComplId>> gc(WITH_HALO);
    // master_printf("Copying the conf from host to device\n");
    // gc=gcH;
    masterPrintf("plaq with new method: %.16lg\n",plaquette(gc));
    // Field<OfComps<Dir,ColorRow,ColorCln,ComplId>,double> gcp(WITH_HALO);


    
//     // auto rrr=(gc*expI).template closeAs<std::decay_t<decltype(gc.createCopy())>>();
//     gc*=expI;
    
//     auto e=(gc*expI).closeAs(gc);
    
    
    
    
//     //gcH=s;
//     master_printf("plaq after phasing: %.16lg\n",plaquette(gc));
//     master_printf("plaq after phasing: %.16lg\n",plaquette(e));
//     // FFT<std::decay_t<decltype(e)>>::fft(e);
//     // master_printf("after fft %lg\n",);
    
// // memcpy(gc.data.storage,gcr._data,sizeof(quad_su3)*locVol);
//     // Field2<OfComps<>,double> t(WITH_HALO);
//     // Field2<OfComps<>,double> u;
//     // t=0.875;
//     // u=shift(t,fw,Dir(0))*shift(t,fw,Dir(0));
//     // printf("u: %lg\n",u(LocLxSite(0)));
    
    
//     // df.updateHalo();
//     // auto fl=df.flatten();
//     // fl(decltype(df)::FlattenedInnerComp(0));
//     // Field2<OfComps<ComplId>,double> ds;
//     // // EoField2<OfComps<SpinRow>,double> df2;
//     // //auto rdf=df.getWritable();
//     // // auto rdf2=df2.getWritable();
//     // df=I;
//     // master_printf("written 0? %lg\n",df(locLxSite(0),spinRow(0),reIm(0)));
//     // master_printf("written 1? %lg\n",df(locLxSite(0),spinRow(0),reIm(1)));
//     // df=df+I;
//     // ds=trace(dag(df(spinRow(0))));
//     // master_printf("written 2? %lg\n",df(locLxSite(0),spinRow(0),reIm(1)));
//     // master_printf("copied -2? %lg\n",ds(locLxSite(0),reIm(1)));
//     // rdf(spinRow(0),locLxSite(0))=1.0;
//     // rdf2(parity(0),spinRow(0),locEoSite(0))=1.0;
//   }
  
// //  crash("");
// //  
// //  auto u=mergedComp<CompsList<Spin,ComplId>>(0);
// //  
// //  /////////////////////////////////////////////////////////////////
// //  
// //  master_printf("allocated in %p\n",conf._data);
// //  {
// //    auto c=conf.getWritable();
// //    master_printf("allocated in %p\n",c._data);
// //    double& e=c[locVol-1][3][2][2][1];
// //    master_printf("end: %p, should be %p\n",&e,c._data+locVol*4*3*3*2-1);
// //  }
// //  
// //  PAR(0,locVol,
// //      CAPTURE(TO_WRITE(conf)),
// //      ivol,
// //      {
// //	su3_put_to_id(conf[ivol][0]);
// //      });
  
//   // start_loc_rnd_gen(235235);
  
//   // spincolor *in=nissa_malloc("in",locVol+bord_vol,spincolor);
//   // spincolor *out=nissa_malloc("out",locVol+bord_vol,spincolor);
//   // spincolor *tmp=nissa_malloc("tmp",locVol+bord_vol,spincolor);
  
//   // /// First test: load a spincolor and unload it
//   // generate_undiluted_source(in,RND_Z2,-1);
  
//   // quda_iface::remap_nissa_to_quda(tmp,in);
//   // quda_iface::remap_quda_to_nissa(out,tmp);
  
//   // master_printf("testing map and unmap, residue: %lg\n",rel_diff_norm(out,in));
  
//   // /// Second test: apply the dirac operator
  
//   // quad_su3 *conf=nissa_malloc("conf",locVol+bord_vol+edge_vol,quad_su3);
//   // spincolor *out_nissa=nissa_malloc("out_nissa",locVol+bord_vol,spincolor);
  
//   // generate_hot_lx_conf(conf);
//   // master_printf("plaq: %lg\n",global_plaquette_lx_conf(conf));
  
//   // vector_reset(in);
//   // in[0][0][0][0]=1.0;
  
//   // const double kappa=0.1325,csw=1.345,mu=0.243;
//   // quda_iface::apply_tmD(out,conf,kappa,csw,mu,in);
  
//   // clover_term_t *Cl=nissa_malloc("Cl",locVol,clover_term_t);
//   // inv_clover_term_t *invCl=nissa_malloc("invCl",locVol,inv_clover_term_t);
//   // clover_term(Cl,csw,conf);
//   // invert_twisted_clover_term(invCl,mu*tau3[0],kappa,Cl);
//   // apply_tmclovQ(out_nissa,conf,kappa,Cl,mu,in);
//   // nissa_free(invCl);
//   // nissa_free(Cl);
  
//   // safe_dirac_prod_spincolor(out_nissa,base_gamma[5],out_nissa);
  
//   // master_printf("comparing\n");
//   // for(int ivol=0;ivol<locVol;ivol++)
//   //   for(int id=0;id<NDIRAC;id++)
//   //     for(int ic=0;ic<NCOL;ic++)
//   // 	for(int ri=0;ri<2;ri++)
//   // 	  {
//   // 	    const double n=out_nissa[ivol][id][ic][ri];
//   // 	    const double q=out[ivol][id][ic][ri];
//   // 	    if(fabs(n)>1e-10 or fabs(q)>1e-10)
//   // 	      master_printf("out,[nissa,quda][ivol=%d,id=%d,ic=%d,ri=%d]: %lg %lg\n",
//   // 			    ivol,id,ic,ri,n,q);
//   // 	  }
//   // master_printf("testing tmD, residue: %lg\n",rel_diff_norm(out,out_nissa));
  
//   // nissa_free(tmp);
//   // nissa_free(in);
//   // nissa_free(out);
//   // nissa_free(out_nissa);
//   // nissa_free(conf);
      
  
  localizer::dealloc();
  fftwFinalize();
      lat.reset();
      delete _lat;
}

int main(int narg,char **arg)
{
  initNissa(narg,arg);
  in_main(narg,arg);
  closeNissa();
  
  return 0;
}
