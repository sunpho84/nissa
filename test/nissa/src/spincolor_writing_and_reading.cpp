#pragma once

int test_spincolor_writing_and_reading()
{
  master_printf("\nGenerating random source of +-1\n");
  
  spincolor *source=nissa_malloc("source",loc_vol,spincolor);
  spincolor *source2=nissa_malloc("source2",loc_vol,spincolor);
  if(nissa_loc_rnd_gen_inited) stop_loc_rnd_gen();
  start_loc_rnd_gen(2342);
  
  generate_undiluted_source(source,RND_Z4,-1);
  write_spincolor("test_wr",source,64);
  read_spincolor(source2,"test_wr");

  int loc_ret=memcmp(source,source2,sizeof(spincolor)*loc_vol);
  
  int ret;
  MPI_Allreduce(&loc_ret,&ret,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  
  if(ret) master_printf("Erorr, read data differs from memory!\n");
  else master_printf("Read the same data present in memory\n");

  master_printf("Removing temporary file \"test_wr\"\n");
  if(rank==0) system("rm -vf test_wr");
  
  nissa_free(source);
  nissa_free(source2);
  
  return !ret;
}
