#include <nissa.hpp>

#include <iostream>

using namespace nissa;

using SpinSpaceColor=TensComps<SpinIdx,LocVolIdx,ColorIdx>;

#ifdef __NVCC__
 #define DEVICE __device__
#else
 #define DEVICE
#endif

// template <typename Fund>
// using SpinColorField=Tens<SpinSpaceColor,Fund>;

// using SpinColorFieldD=SpinColorField<double>;

/*
#define TEST(NAME,...)							\
   double& NAME(SpinColorFieldD& tensor,SpinIdx spin,ColorIdx col,LocVolIdx space) \
   {									\
     asm("#here " #NAME "  access");					\
     return __VA_ARGS__;							\
   }
   */

// LocVolIdx vol;

// struct _Spin : public TensCompSize<4>
// {
// };

// TEST(seq_fun,bindComp(tensor,col,spin,space).eval())

// TEST(bra_fun,tensor[col][spin][space])

// TEST(csv_fun,tensor(col,spin,space));

// TEST(svc_fun,tensor(spin,space,col));

// TEST(hyp_fun,tensor(col)(spin)(space));

// TEST(triv_fun,tensor.trivialAccess(col+3*(space+vol*spin)));

// template <typename T>
// void test_if_ref(T&& t)
// {
//   const bool b=is_const_lvalue_reference_v<T>;
//   std::cout<<"ECCO "<<b<<std::endl;
// }

// template <typename TOut,
// 	  typename TIn1,
// 	  typename TIn2>
// void unsafe_su3_prod(TOut&& out,const TIn1& in1,const TIn2& in2)
// {
//   LocVolIdx s(0);
  
//   for(RwColorIdx i1{0};i1<3;i1++)
//     for(ClColorIdx k2{0};k2<3;k2++)
//       {
//   	out(i1,k2,s)=0.0;
//   	for(ClColorIdx i2(0);i2<3;i2++)
//   	  out(i1,k2,s)+=in1(i1,i2,s)*in2(i2.transp(),k2,s);
//       }
// }

using G=  Field<LocVolIdx,SpinColorComplComps,double,FieldLayout::GPU>;
using C=  Field<LocVolIdx,SpinColorComplComps,double,FieldLayout::CPU>;

void maybe(G& testG,C& testC)
{
  asm("#here starts maybe");
  G::spaceTimeParallelLoop([testC=testC.data.getRef(),testG=testG.data.getRef()] DEVICE
  					 (const int _ivol) mutable
  					 {
					   LocVolIdx ivol{_ivol};
  					   for(SpinIdx id{0};id<4;id++)
  					     for(ColorIdx ic{0};ic<3;ic++)
  					       for(ComplIdx ri{0};ri<2;ri++)
  						 {
  						   // double* offC=testC.data.data.data;
  						   // double* offG=testG.data.data.data;
  						   double& tC=testC(ivol,id,ic,ri);
  						   double& tG=testG(ivol,id,ic,ri);
  						   // int nC=(int)((size_t)&tC-(size_t)offC)/sizeof(double);
  						   // int nG=(int)((size_t)&tG-(size_t)offG)/sizeof(double);
  						   // printf("ivol %d loc_vol %ld %d %d\n",(int)ivol,loc_vol,nC,nG);
  						   tG=tC;
  						 }
  					 });
  asm("#here ends maybe");
}

void trivial(G& testG,C& testC)
{
  asm("#here starts again");
  double* tC=testC.data.data.data.data;
  double* tG=testG.data.data.data.data;;
  
  for(Size ivol=0;ivol<loc_vol;ivol++)
    for(Size id=0;id<4;id++)
      for(Size ic=0;ic<3;ic++)
	for(Size ri=0;ri<2;ri++)
	  tG[ivol+loc_vol*(ri+2*(ic+3*id))]=tC[ri+2*(ic+3*(id+4*ivol))];
  
  asm("#here ends again");
  master_printf("%lg\n",testG(LocVolIdx{3},SpinIdx{1},ColorIdx{2},ComplIdx{0}));
}


void in_main(int narg,char **arg)
{
  init_grid(8,4);
  
  Field<LocVolIdx,SpinColorComplComps,double,FieldLayout::GPU,TensStorageLocation::ON_GPU> testG_onG(HaloKind::NO_HALO);
  
  G::spaceTimeParallelLoop([testG_onG=testG_onG.data.getRef()] DEVICE
  					 (const int _ivol) mutable
  					 {
					   LocVolIdx ivol{_ivol};
  					   for(SpinIdx id{0};id<4;id++)
  					     for(ColorIdx ic{0};ic<3;ic++)
  					       for(ComplIdx ri{0};ri<2;ri++)
  						 {
  						   double& tG=testG_onG(ivol,id,ic,ri);
  						   tG=M_PI;
  						 }
  					 });
  
  //temporary hack
  Field<LocVolIdx,SpinColorComplComps,double,FieldLayout::GPU,TensStorageLocation::ON_CPU> testG_onC(HaloKind::NO_HALO);
  decript_cuda_error(cudaMemcpy(testG_onC.data.data.data.data,testG_onG.data.data.data.data,sizeof(spincolor)/sizeof(double)*loc_vol,cudaMemcpyDeviceToHost),"Copying");
  //gpu_to_cpu_copy(testG_onC,testG_onG);
  
  for(LocVolIdx ivol{0};ivol<loc_vol;ivol++)
    for(SpinIdx id{0};id<4;id++)
      for(ColorIdx ic{0};ic<3;ic++)
	for(ComplIdx ri{0};ri<2;ri++)
	  master_printf("\n",testG_onC(ivol,id,ic,ri));
  
  // trivial(testG,testC);
}

int main(int narg,char **arg)
{
//   using std::cout;
//   using std::cin;
//   using std::endl;
  
//   cout<<"Please enter volume: ";
//   cin>>vol;
  
//   /// Spinspacecolor
//   Tens<SpinSpaceColor,double> tensor(vol);
  
//   /// Fill the spincolor with flattened index
//   for(SpinIdx s(0);s<4;s++)
//     for(LocVolIdx v(0);v<vol;v++)
//       for(ColorIdx c(0);c<3;c++)
//   	tensor(s,v,c)=c+3*(v+vol*s);
  
//   // Read components to access
//   cout<<"Please enter spin, space and color index to printout: ";
  
//   /// Spin component
//   SpinIdx spin;
//   cin>>spin;
  
//   /// Space component
//   LocVolIdx space;
//   cin>>space;
  
//   /// Color component
//   ColorIdx col;
//   cin>>col;
  
//   double& hyp=hyp_fun(tensor,spin,col,space);
  
//   /// Spin,space,color access
//   double& svc=svc_fun(tensor,spin,col,space);
  
//   /// Color,spin,space access
//   double& csv=csv_fun(tensor,spin,col,space);
  
//   /// Color,spin,space access
//   double& seq=seq_fun(tensor,spin,col,space);
  
//   /// Color,spin,space access with []
//   double& bra=bra_fun(tensor,spin,col,space);
  
//   /// Trivial spin access
//   double& t=triv_fun(tensor,spin,col,space);
  
//   using SU3FieldComps=TensComps<RwColorIdx,ClColorIdx,LocVolIdx>;
  
//   using SU3Field=Tens<SU3FieldComps,double>;
  
//   SU3Field conf1(vol),conf2(vol),conf3(vol);
  
//   using ComplComps=TensComps<ComplIdx>;
//   // Binder<const Tens<std::tuple<TensorCompIdx<_Color, ROW, 0>,
//   // 			       TensorCompIdx<_Color, COL, 0>,
//   // 			       TensorCompIdx<_Space, ANY, 0> >,
//   // 		    double>&,
//   // 	 const Tens<std::tuple<TensorCompIdx<_Color, ROW, 0>,
//   // 			       TensorCompIdx<_Color, COL, 0>,
//   // 			       TensorCompIdx<_Space, ANY, 0> >,
//   // 		    double>&,
//   // 	 TensorCompIdx<_Color, COL, 0>,
//   // 	 TensorCompIdx<_Color, ROW, 0>,
//   // 	 TensorCompIdx<_Space, ANY, 0> >
//   //const double a=0.0;
//   //remove_const_if_ref(a)=0.0;
//     Tens<ComplComps,double> test;
//   test.trivialAccess(0)=0.0;
  
//     for(LocVolIdx v(0);v<vol;v++)
//       for(RwColorIdx c1(0);c1<3;c1++)
// 	for(ClColorIdx c2(0);c2<3;c2++)
// 	  conf1(space,c1,c2)=conf2(space,c1,c2)=conf3(space,c1,c2)=0.0;
  
//   unsafe_su3_prod(conf1,conf2,conf3);
  
//   conf1(ClColorIdx{0},RwColorIdx{1},LocVolIdx{0})=1.0;
//   cout<<"Transp: "<<transpose(conf1(LocVolIdx{0})(RwColorIdx{1}))(RwColorIdx{0})<<endl;
  
//   cout<<"Test:";
//   cout<<" "<<conf1[space][col][col.transp()];
//   cout<<" "<<svc;
//   cout<<" "<<csv;
//   cout<<" "<<t;
//   cout<<" "<<seq;
//   cout<<" "<<hyp;
//   cout<<" "<<bra;
//   cout<<" expected: "<<col+3*(space+vol*spin)<<endl;
  
//   cout<<SpinColorFieldD::stackAllocated<<endl;
//   cout<<SU3Field::stackAllocated<<endl;
  
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
