#pragma once

int test_random_source_generation()
{
  master_printf("\nGenerating random source of +-1\n");
  
  spincolor *source=appretto_malloc("source",loc_vol,spincolor);
  start_loc_rnd_gen(2342);
  
  generate_undiluted_source(source,RND_Z4,-1);
  checksum check_save={1857847612,766990187};
    
  checksum check_comp;
  checksum_compute_appretto_data(check_comp,source,sizeof(spincolor));
  
  master_printf("checksum computed: %u %u\n",check_comp[0],check_comp[1]);
  master_printf("checksum saved:    %u %u\n",check_save[0],check_save[1]);
  
  appretto_free(source);
  
  return 1;
}
