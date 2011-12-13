#pragma once

extern "C" void test_dirac_operator(spincolord *cpu_result,quad_su3d *cpu_uncomp_conf,spincolord *cpu_source,double mass,double kappa)
{
  //Compress the configuration on cpu and send it to the gpu
  su3c *gpu_comp_conf_high=gpu_allocate_su3c(2*4*(loc_vol+loc_bord),"gpu_comp_conf");
  su3c *gpu_comp_conf_low=gpu_comp_conf_high+4*2*(loc_vol+loc_bord);
  
  gaugeconf_compress_and_move(gpu_comp_conf_high,gpu_comp_conf_low,cpu_uncomp_conf,loc_vol+loc_bord);
  
  //Compress the spincolor on the cpu and send it to the gpu
  
  //free the gpu compressed conf and cpu decompressed conf
  gpu_free(gpu_comp_conf_high,"gpu_comp_conf");
}
