#pragma once

void back()
{
  if(rank==0) printf("Backtracing...\n");
  void *callstack[128];
  int i,frames=backtrace(callstack,128);
  char **strs=backtrace_symbols(callstack,frames);
  if(rank==0) for(i=0;i<frames;i++) printf("%s\n",strs[i]);

  free(strs);               
}

void terminate_sigsegv(int par)
{
  print_all_appretto_vect_content();
  back();
  crash("Error detected");
}
