#pragma once

extern "C" void test_dirac_operator(spincolord *cpu_result,quad_su3d *cpu_uncomp_conf,spincolord *cpu_source,double mass,double kappa)
{
  //Compress the configuration on cpu and send it to the gpu
  quad_su3c *gpu_comp_conf=gpu_allocate_quad_su3c(loc_vol+loc_bord,"gpu_comp_conf");
  gaugeconf_compress_and_move(gpu_comp_conf,conf,loc_vol+loc_bord);
  
  //Compress the spincolor on the cpu and send it to the gpu
  
  //free the gpu compressed conf and cpu decompressed conf
  gpu_free(gpu_comp_conf,"gpu_comp_conf");
}
