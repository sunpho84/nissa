#ifndef _DEBUG_H
#define _DEBUG_H
double take_time();
void internal_crash(int line,const char *file,const char *templ,...);
void internal_crash_printing_error(int line,const char *file,int err_code,const char *templ,...);
void internal_decript_MPI_error(int line,const char *file,int rc,const char *templ,...);
void print_backtrace_list();
void terminate_sigsegv(int);
#endif
