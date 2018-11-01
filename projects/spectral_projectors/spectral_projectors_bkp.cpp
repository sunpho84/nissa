#include "nissa.hpp"

using namespace nissa;

void init_test(int argc, char** argv)
{
    init_nissa(argc,argv);
  
  //init the grid
  init_grid(8,4);
}

void test(int passed,const char *test_name)
{
  if(!passed) master_printf("\n%s test not passed!\n\n",test_name);
  else master_printf("\n%s test passed\n\n",test_name);
  
  master_printf("################ %s test finished ###############\n",test_name);
}


int test_plaquette_computation()
{
  const double exp_plaq=0.6350055722288719;
  const double tolerance=1.e-14;
  
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //read conf, compute plaquette, print it
  read_ildg_gauge_conf(conf,"../../data/L4T8conf");
  
  double plaq=global_plaquette_lx_conf(conf);
  master_printf("Loaded plaquette: %16.16g, expected: %16.16g\n",plaq,exp_plaq);
  
  //check precision
  double rel_diff=(plaq-exp_plaq)/exp_plaq;
  int test_passed=fabs(rel_diff)<=tolerance;  
  master_printf("Relative difference: %16.16g\n",rel_diff);
  master_printf("Tolerance: %16.16g\n",tolerance);
  
  nissa_free(conf);
  
  return test_passed;
}

//  //find the neig eigenvalues closest to the target
//    template <class Fmat,class Filler>
//    void eigenvalues_of_hermatr_find_autarchic(complex **eig_vec,complex *eig_val,int neig,bool min_max,
//  					     const int mat_size,const int mat_size_to_allocate,const Fmat &imp_mat,
//  					     const double target_precision,const int niter_max,
//  					     const Filler &filler)

//wrap the generation of the test vector into an object that can be passed to the eigenfinder
const auto filler=[](complex *a)
{
  generate_undiluted_source((spincolor*)a,RND_GAUSS,ALL_TIMES);
};


int test_eigenvalue_computation()
{
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //read conf, compute first 10 eigenvalues and eigenvectors of D^+D, printing eigenvalues
  read_ildg_gauge_conf(conf,"../../data/L4T8conf");

  //std::vector<double> eigs = compute_DD_eigenvalues(conf);  
  eigenvalues_of_hermatr_find((complex**)eigvec,Q2_eig_val,neig,min_max,mat_size,mat_size_to_allocate,imp_mat,eig_precision,niter_max,filler);



//  double plaq=global_plaquette_lx_conf(conf);
//  master_printf("Loaded plaquette: %16.16g, expected: %16.16g\n",plaq,exp_plaq);
//  
//  //check precision
//  double rel_diff=(plaq-exp_plaq)/exp_plaq;
//  int test_passed=fabs(rel_diff)<=tolerance;  
//  master_printf("Relative difference: %16.16g\n",rel_diff);
//  master_printf("Tolerance: %16.16g\n",tolerance);
  
  nissa_free(conf);
  
  return 0;
}


int main(int argc, char** argv){
  init_test(argc, argv);
  
  test(test_plaquette_computation(),"Plaquette test");
  

  test(test_eigenvalue_computation(), "Eigenvalue test");


  close_nissa();
  
  return 0;
}
