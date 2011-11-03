#pragma once

int test_write_and_read()
{
  master_printf("\nGenerating random source of +-1\n");
  
  spincolor *source=appretto_malloc("source",loc_vol,spincolor);
  spincolor *source2=appretto_malloc("source2",loc_vol,spincolor);
  start_loc_rnd_gen(2342);
  
  generate_undiluted_source(source,RND_Z4,-1);
  write_spincolor("test_wr",source,64);
  read_spincolor(source2,"test_wr");

  int loc_ret=memcmp(source,source2,sizeof(spincolor)*loc_vol);
  
  int ret;
  MPI_Allreduce(&loc_ret,&ret,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  
  if(ret) master_printf("Erorr, read data differs from memory!\n");
  else master_printf("Read the same data than in memory\n");
  
  if(rank==0) system("rm -vf test_wr");
  
  appretto_free(source);
  appretto_free(source2);
  
  return ret;
}
