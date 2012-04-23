void internal_crash(int line,const char *file,const char *templ,...);
void internal_decript_MPI_error(int line,const char *file,int rc,const char *templ,...);
void print_backtrace_list();
void terminate_sigsegv(int par);
