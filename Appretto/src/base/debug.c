#pragma once

void print_backtrace_list()
{
  void *callstack[128];
  int frames=backtrace(callstack,128);
  char **strs=backtrace_symbols(callstack,frames);
  
  master_printf("Backtracing...\n");
  for(int i=0;i<frames;i++) master_printf("%s\n",strs[i]);
  
  free(strs);               
}

void terminate_sigsegv(int par)
{
  print_all_appretto_vect_content();
  print_backtrace_list();
  crash("Error detected");
}
