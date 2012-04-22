#include "src/common.cpp"
#include "src/plaquette_computation.cpp"
#include "src/random_source_generation.cpp"
#include "src/spincolor_writing_and_reading.cpp"
#include "src/Q2tm_cg_inversion.cpp"

int main()
{
  init_test();
  
  test(test_plaquette_computation(),"plaquete computation");
  test(test_random_source_generation(),"random source generation");
  test(test_spincolor_writing_and_reading(),"writing and reading random spinor");
  test(test_Q2tm_inversion(),"Q2tm inversion");
  
  close_nissa();
  
  return 0;
}
