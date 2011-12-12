#pragma once

int test_random_source_generation()
{
  master_printf("\nGenerating random source of +-1\n");
  
  spincolor *source=nissa_malloc("source",loc_vol,spincolor);
  if(nissa_loc_rnd_gen_inited) stop_loc_rnd_gen();
  start_loc_rnd_gen(2342);
  
  generate_undiluted_source(source,RND_Z4,-1);
  checksum check_save={0x5a0f8506,0x576e64ff};
  
  checksum check_comp;
  checksum_compute_nissa_data(check_comp,source,sizeof(spincolor));
  
  master_printf("Checksum computed: %#010x %#010x\n",check_comp[0],check_comp[1]);
  master_printf("Checksum saved:    %#010x %#010x\n",check_save[0],check_save[1]);
  
  nissa_free(source);
  
  return ((check_comp[0]==check_save[0])&&(check_comp[1]==check_save[1]));
}
